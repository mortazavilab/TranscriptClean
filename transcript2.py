# This file contains the Transcript2 class and associated methods for the 
# TranscriptClean program

from spliceJunction import SpliceJunction
from intronBound import IntronBound
import pyfasta
import pybedtools
import re
import itertools
import string

class Transcript2:
    def __init__(self, sam, genome, spliceAnnot):
        samFields = sam.strip().split('\t')

        # Keep track of changes for log file
        #self.log = []
        #self.transcriptErrors = []

        # These eleven attributes are initialized directly from the input 
        # SAM entry and are mandatory 
        self.QNAME = samFields[0]
        self.FLAG = samFields[1]
        self.CHROM = samFields[2]
        self.POS = int(samFields[3])
        self.MAPQ = samFields[4]
        self.CIGAR = samFields[5]
        self.RNEXT = samFields[6]
        self.PNEXT = samFields[7]
        self.TLEN = samFields[8]
        self.SEQ = samFields[9]
        self.QUAL = "*"

        self.logInfo = [self.QNAME, "primary", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] 
        # Check if read is unmapped (0), uniquely mapped (1), multimapped (2)
        self.mapping = 1
        if self.CHROM == "*" or int(self.FLAG) == 4: 
            self.mapping = 0
            self.logInfo = [self.QNAME, "unmapped", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
        elif int(self.FLAG) > 16:
            self.mapping = 2
            self.logInfo = [self.QNAME, "non-primary", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]

        # If the sam entry contains additional optional fields, process them 
        self.NM = ""
        self.MD = ""
        self.jM = ""
        self.jI = ""        
        otherFields = []

        for field in samFields[11:len(samFields)]:
            if field.startswith("NM"): self.NM = field
            elif field.startswith("MD"): self.MD = field
            elif field.startswith("jM"): self.jM = field 
            elif field.startswith("jI"): self.jI = field
            else: otherFields.append(field)

        # Compute NM and MD tags here.
        if self.mapping == 1:
            self.NM, self.MD = self.getNMandMDFlags(genome)

        # If the NM and MD tags are None, it means there was a reference genome
        # problem somewhere in the read. Consider such reads unmapped.
        if self.NM == None and self.MD == None:
            self.mapping = 0
            self.logInfo = [self.QNAME, "unmapped", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]

        # If the jM and jI fields are missing, compute them here.
        if (self.jM == self.jI == "") and self.mapping == 1:
            self.jM, self.jI = self.getjMandjITags(genome, spliceAnnot)

        self.otherFields = "\t".join(otherFields)        

        # These attributes are set by parsing the inputs
        self.spliceJunctions = []
        self.isCanonical = True
        self.strand = "+"        
        if int(self.FLAG) == 16 or int(self.FLAG) == 272: 
            self.strand = "-"

        # Only run this section if there are splice junctions
        if self.jM != "" and "-1" not in self.jM:
            # Create an object for each splice junction
            self.spliceJunctions = self.parseSpliceJunctions(genome)            

    def addCorrectedDeletion(self):
        """ Add a corrected deletion to log"""
        self.logInfo[2] += 1
        return    

    def addUncorrectedDeletion(self):
        """ Add an uncorrected deletion to log"""
        self.logInfo[3] += 1
        return

    def addVariantDeletion(self):
        """ Add a variant deletion to log"""
        self.logInfo[4] += 1
        return

    def addCorrectedInsertion(self):
        """ Add a corrected insertion to log"""
        self.logInfo[5] += 1
        return    

    def addUncorrectedInsertion(self):
        """ Add an uncorrected insertion to log"""
        self.logInfo[6] += 1
        return

    def addVariantInsertion(self):
        """ Add a variant insertion to log"""
        self.logInfo[7] += 1
        return

    def addCorrectedMismatch(self):
        """ Add a corrected mismatch to log"""
        self.logInfo[8] += 1
        return

    def addVariantMismatch(self):
        """ Add a variant mismatch to log"""
        self.logInfo[9] += 1
        return

    def addCorrected_NC_SJ(self):
        """ Add a corrected noncanonical SJ to log"""
        self.logInfo[10] += 1
        return

    def addUncorrected_NC_SJ(self):
        """ Add an uncorrected noncanonical SJ to log"""
        self.logInfo[11] += 1
        return
    
    def recheckCanonical(self):
        """ Check each splice junction. If one or more junctions are
            noncanonical, then so is the transcript. """
        for jn in self.spliceJunctions:
            if jn.isCanonical == False:
                self.isCanonical = False
                return False
        self.isCanonical = True
        return True

    def compute_transcript_end(self):
        """ Given the start position and CIGAR string of a mapped SAM transcript,
            compute the end position in the reference genome.

            Args:
                start: The start position of the transcript with respect to the
                forward strand

                cigar: SAM CIGAR string describing match operations to the reference
                genome

            Returns:
                end position of the transcript.
        """
        end = self.POS

        ops, counts = self.splitCIGAR()
        for op,ct in zip(ops, counts):
            if op in ["M", "N", "D"]:
                end += ct

        return end - 1

    def splitCIGAR(self):
        """ Takes CIGAR string from SAM and splits it into two lists: 
            one with capital letters (match operators), and one with 
            the number of bases that each operation applies to. """

        alignTypes = re.sub('[0-9]', " ", self.CIGAR).split()
        counts = re.sub('[A-Z]', " ", self.CIGAR).split()
        counts = [int(i) for i in counts]

        return alignTypes, counts

    def splitMD(self):
        """ Takes MD tag and splits into two lists: 
            one with capital letters (match operators), and one with 
            the number of bases that each operation applies to. """

        MD = str(self.MD).split(":")[2]
        operations = []

        # Split MD string where type changes. 
        # Digits are separated from base changes. 
        # Deletions (with ^) are captured together.
        counts = ["".join(x) for _, x in itertools.groupby(MD, key=str.isdigit)]
        
        # Get operations
        for i in range(0,len(counts)):
            curr = counts[i]
            try:
                counts[i] = int(curr)
                operations.append("M")
            except ValueError:
                # Handle deletion 
                if curr.startswith("^"): 
                    operations.append("D")
                    counts[i] = len(counts[i]) - 1
                else: 
                    operations.append("X")
                    counts[i] = len(counts[i])
       
        return operations, counts

    def mergeMDwithCIGAR(self):
        """ Takes the MD and CIGAR strings, and combines them into a unified 
            structure that encodes all possible operations w.r.t the reference:
            match, mismatch, deletion, insertion, hard clipping, 
            and soft clipping. """

        mergeCounts = []
        mergeOperations = []
    
        cigarOperation, cigarCount = self.splitCIGAR()
        mdOperation, mdCount = self.splitMD() 

        mdIndex = 0
        cigarIndex = 0

        while mdIndex < len(mdOperation) or cigarIndex < len(cigarOperation):

            # If the current CIGAR operation is S, H, N, or I, add that to the
            # output. The MD tag doesn't have these
            if cigarOperation[cigarIndex] in ("H", "S", "I", "N"):
                    mergeOperations.append(cigarOperation[cigarIndex])
                    mergeCounts.append(cigarCount[cigarIndex])
                    cigarIndex += 1

            # Otherwise, we need to compare the current CIGAR and MD operations.
            # Select the "shorter" operation and add it to the results. 
            # Subtract away the same number of bases from the competing entry.
            else:
                if cigarCount[cigarIndex] < mdCount[mdIndex]:
                # If the CIGAR string lists fewer matched bases than MD, 
                # it means the CIGAR has had an insertion not listed in MD
                    mdCount[mdIndex] = mdCount[mdIndex] - cigarCount[cigarIndex]
                    mergeOperations.append(cigarOperation[cigarIndex])
                    mergeCounts.append(cigarCount[cigarIndex])
                    cigarIndex += 1

                elif cigarCount[cigarIndex] > mdCount[mdIndex]:
                # If the CIGAR string lists more matched bases than MD, 
                # it means that MD has a mismatch not listed in CIGAR
                    cigarCount[cigarIndex] = cigarCount[cigarIndex] - mdCount[mdIndex]
                    mergeOperations.append(mdOperation[mdIndex])
                    mergeCounts.append(mdCount[mdIndex])
                    mdIndex += 1
                    
                # For cases where both MD and CIGAR specify the same match type, 
                # add to the result and advance to next position in lists
                else: 
                    mergeOperations.append(mdOperation[mdIndex])
                    mergeCounts.append(mdCount[mdIndex])
                    mdIndex += 1
                    cigarIndex += 1

        return mergeOperations, mergeCounts


    def parseSpliceJunctions(self, genome):
        """ Takes the splice junction information from the SAM input and 
            creates a SpliceJunction object for each junction."""

        spliceJns = ((self.jM).split(":")[-1]).split(",")[1:]
        intronBounds = ((self.jI).split(":")[-1]).split(",")[1:]

        count = 0
        jnNum = 0
        jnObjects = [] 
        for entry in spliceJns:
            start = int(intronBounds[count])
            end = int(intronBounds[count + 1])
            sj = SpliceJunction(self.QNAME, jnNum, self.CHROM, start, end, \
                 self.strand, entry, genome)
            jnObjects.append(sj)

            # Check if junction is canonical or not. 
            if sj.isCanonical == False: self.isCanonical = False
            count += 2
            jnNum += 1
        
        return jnObjects

    def printableSAM(self, genome, spliceAnnot):
        """ Returns a SAM-formatted string representation of the transcript"""
        if len(self.spliceJunctions) > 0:
            self.jI = "jI:B:i," + ",".join(str(i.pos) for i in self.getAllIntronBounds())
            self.jM = "jM:B:c," + ",".join(str(i) for i in self.getAllSJMotifs(genome, spliceAnnot))
        self.NM, self.MD = self.getNMandMDFlags(genome)        
        fields = [ self.QNAME, self.FLAG, self.CHROM, self.POS, self.MAPQ, self.CIGAR, \
                   self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, self.QUAL, self.otherFields, \
                   self.NM, self.MD, self.jM, self.jI ]
        fields = filter(None, fields)
        return "\t".join([str(x) for x in fields]).strip()

    def printableFa(self):
        """ Returns a fasta-formatted string representation of the transcript """
        fastaID = ">" + self.QNAME
        strand = self.strand
        seq = str(self.SEQ)

        if strand == "-": # Need to reverse-complement the sequence
            seq = reverseComplement(seq)
   
        # Split seq into 80-character segments 
        fastaSeq = [seq[i:i+80] for i in range(0, len(seq), 80)]
        return fastaID + "\n" + "\n".join(fastaSeq)

    def getAllIntronBounds(self):
        """ Return all intron bound objects belonging to this transcript """

        result = []
        for jn in self.spliceJunctions:
            b = jn.bounds
            result.append(b[0])
            result.append(b[1])
        return result
   
    def getAllSJMotifs(self, genome, spliceAnnot):
        """ Return all splice junction motifs translated into their numeric 
            STAR codes"""
        result = []
        for jn in self.spliceJunctions:
            #SpliceJunction.recheckJnStr(jn, genome, spliceAnnot)
            result.append(jn.jnStr)
        return result
 
    def getNMandMDFlags(self, genome):
        """ This function uses the transcript sequence, its CIGAR string, 
            and the reference genome to create NM and MD sam flags."""
        NM = 0
        MD = "MD:Z:"
        MVal = 0
        seqPos = 0
        genomePos = self.POS

        if self.mapping != 1:
            return "",""

        operations, counts = self.splitCIGAR()

        tot = 0
        for op, ct in zip(operations, counts):
            if op in ["M", "I", "S"]:
                tot += ct
        
        for op, ct in zip(operations, counts):
            if op == "M":
                for i in range(0,ct):
                    currBase = self.SEQ[seqPos]
                    refBase = genome.sequence({'chr': self.CHROM, 'start': genomePos, 'stop': genomePos}, one_based=True) 
                    if refBase == "":
                        return None, None
                   
                    # In the event of a mismatch
                    if currBase.upper() != refBase.upper():
                        # End any match we have going and add the mismatch
                        MD = MD + str(MVal)  
                        MVal = 0 
                        MD = MD + str(refBase) 
                        NM += 1
                    # Bases match
                    else:
                        MVal += 1
                    # Either way, advance forwards in the sequence and genome
                    seqPos += 1
                    genomePos += 1
            if op == "D":
                # End any match we have going and add the missing reference bases
                MD = MD + str(MVal)  
                MVal = 0
                refBases = genome.sequence({'chr': self.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                if refBase == "":
                    return None, None
                MD = MD + "^" + str(refBases)
                NM += ct
                genomePos += ct
            # For insertions and soft clips, we move on without adding to the MD
            if op in ["I", "S"]:
                seqPos += ct
                if op == "I": NM += ct
            if op in ["N", "H"]:
                genomePos += ct
            
        if MVal > 0: MD = MD + str(MVal) 
        return "NM:i:" + str(NM), MD

    def getjMandjITags(self, genome, spliceAnnot):
        """ If the input sam file doesn't have the custom STARlong-derived jM 
            and jI tags, we need to compute them. This is done by stepping 
            through the CIGAR string and sequence. When an intron (N) is 
            encountered, we check the first two bases and last two bases of 
            the intron in the genome sequence to detemine whether they are 
            canonical. We also record the start and end position of the intron. """

        if self.mapping != 1:
            return "",""
       
        seq = self.SEQ
        operations, counts = self.splitCIGAR()

        jM = ["jM:B:c"] 
        jI = ["jI:B:i"]
        
        genomePos = self.POS

        # Iterate over operations
        for op,ct in zip(operations, counts):
            if op == "N":
                # This is an intron
                intronStart = genomePos
                startBases = genome.sequence({'chr': self.CHROM, 'start': genomePos, 'stop': genomePos + 1}, one_based=True)
                intronEnd = genomePos + ct - 1
                endBases = genome.sequence({'chr': self.CHROM, 'start': intronEnd - 1, 'stop': intronEnd}, one_based=True)
             
                # Check if junction is annotated
                if (self.CHROM + "_" + str(intronStart)) in spliceAnnot and (self.CHROM + "_" + str(intronEnd)) in spliceAnnot:
                    motifCode = 20 + getSJMotifCode(startBases, endBases)
                else: 
                    motifCode = getSJMotifCode(startBases, endBases)
 
                jM.append(str(motifCode))
                jI.append(str(intronStart))
                jI.append(str(intronEnd))

            if op not in ["S", "I"]:
                 genomePos += ct

        # If the transcript has no introns, we need to add -1 to the tags
        if len(jM) == len(jI) == 1:
            jM.append("-1")
            jI.append("-1")

        jMstr = ",".join(jM)
        jIstr = ",".join(jI)

        return jMstr, jIstr

    def base_wise_CIGAR(self):
        """ Create an extended version of the CIGAR string with an operation
            per base.
            Example: 3M10N2D5M becomes
                     MMMNNNNNNNNNNDDMMMMM

        """
        base_wise_str = ""
        ops, counts = self.splitCIGAR()        
        for op,ct in zip(ops, counts):
            base_wise_str += op*ct
        return base_wise_str
    
    def fetch_region_sequence(self, chromosome, start, end):
        """ Walks the SAM sequence to return the bases in the specified region.
            Returns None if the sequence is not available, i.e. because the 
            region does not overlap with the transcript, or because it is in
            an intron. Supplied coordinates should be 1-based. """
         
        # Check whether the supplied region is located even remotely near the
        # transcript
        if chromosome != self.CHROM:
            return None
        if not(start >= self.POS and end <= self.compute_transcript_end()):
            return None
        
        # Walk transcript sequence using CIGAR string
        positions = range(start, end + 1)
        seq = self.SEQ
        seq_pos = 0
        bases = ""

        genome_pos = self.POS
        while genome_pos <= end:
            for op in self.base_wise_CIGAR():
                if op == "M":
                    if (genome_pos in positions):
                        bases += seq[seq_pos]
                    genome_pos += 1
                    seq_pos += 1
                # Advance in genome sequence but not in transcript sequence
                if op in ["D", "N", "H"]:
                    if (genome_pos in positions) and op == "D":
                        bases += "-"
                    genome_pos += 1
                # Advance in transcript sequence but not genome sequence
                if op in ["S", "I"]:
                    if (genome_pos in positions) and op == "I":
                        bases += seq[seq_pos]
                    seq_pos += 1       
 
        if bases == "":
            return None

        # If the transcript is on the reverse strand, reverse-complement the
        # sequence before returning it
        if self.strand == "-":    
            bases = reverseComplement(bases)
        return bases
        

def getSJMotifCode(startBases, endBases):
    """ Determines which STAR-style splice junction code applies to a splice motif """       

    motif = (startBases + endBases).upper()

    if motif == "GTAG":
        return 1
    elif motif == "CTAC":
        return 2
    elif motif == "GCAG":
        return 3
    elif motif == "CTGC":
        return 4
    elif motif == "ATAC":
        return 5
    elif motif == "GTAT":
        return 6
    else:
        return 0
        
def reverseComplement(seq):
    """ Returns the reverse complement of a DNA sequence, 
        retaining the case of each letter"""
    complement = ""

    for base in seq:
        if base == "A": complement += "T"
        elif base == "T": complement += "A" 
        elif base == "G": complement += "C"
        elif base == "C": complement += "G"
        elif base == "N": complement += "N"
        elif base == "a": complement += "t"
        elif base == "t": complement += "a" 
        elif base == "g": complement += "c"
        elif base == "c": complement += "g"
	elif base == "n": complement += "n"
        elif base == "*": complement += "*"
        else:
            complement += base
            #print "Warning: reverse complement function encountered unknown base " + "'" + base + "'"

    reverseComplement = complement[::-1]

    return reverseComplement

