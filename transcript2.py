# This file contains classes for the TranscriptClean program

from spliceJunction import SpliceJunction
from intronBound import IntronBound
import pyfasta
import pybedtools
import re
import itertools

class Transcript2:

    def __init__(self, sam, genome):
        samFields = sam.strip().split('\t')

        # These eleven attributes are initialized directly from the input SAM entry and are mandatory 
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
        self.QUAL = samFields[10]

        # If the sam entry contains additional optional fields, process them here
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

        # If the NM and MD tags were missing, compute them here.
        if self.MD == "":
            self.NM, self.MD = self.getNMandMDFlags(genome)

        self.otherFields = "\t".join(otherFields)        

        # These attributes are set by parsing the inputs
        self.spliceJunctions = []
        self.isCanonical = True
        self.strand = "+"        
        if self.FLAG == 16: self.strand = "-"

        # Only run this section if there are splice junctions
        if self.jM != "" and "-1" not in self.jM:
            # Create an object for each splice junction
            self.spliceJunctions = self.parseSpliceJunctions(genome)            
        
    def recheckCanonical(self):
        for jn in self.spliceJunctions:
            if jn.isCanonical == False:
                self.isCanonical = False
                return False
        self.isCanonical = True
        return True


    def splitCIGAR(self):
        # Takes CIGAR string from SAM and splits it into two lists: one with capital letters (match operators), and one with the number of bases

        alignTypes = re.sub('[0-9]', " ", self.CIGAR).split()
        counts = re.sub('[A-Z]', " ", self.CIGAR).split()
        counts = [int(i) for i in counts]

        return alignTypes, counts

    def splitMD(self):
        # Takes MD tag and splits into individual operations

        MD = self.MD.split(":")[2]
        operations = []

        # Split MD string where type changes. Digits are separated from base changes. Deletions (with ^) are captured together.
        counts = ["".join(x) for _, x in itertools.groupby(MD, key=str.isdigit)]

        # Get operations
        for i in range(0,len(counts)):
            curr = counts[i]
            try:
                counts[i] = int(curr)
                operations.append("M")
            except ValueError:
                #Handle the exception
                if curr.startswith("^"): 
                    operations.append("D")
                    counts[i] = len(counts[i]) - 1
                else: 
                    operations.append("X")
                    counts[i] = len(counts[i])

        return operations, counts

    def mergeMDwithCIGAR(self):
        # This function takes the MD and CIGAR strings, and combines them into a unified structure that encodes all possible operations w.r.t the reference: match, mismatch, deletion, insertion, hard clipping, and soft clipping.

        mergeCts = []
        mergeOps = []
    
        cOp, cCt = self.splitCIGAR()
        mOp, mCt = self.splitMD()  

        cCt = [31, 1, 17, 1, 37]
        cOp = ["M", "I", "M", "D", "M"]

        mCt = [6,1,4,1,20,1,1,1,5,1,5,1,1,1,3,1,15,1,1,1,15]
        mOp = ["M", "X", "M", "X", "M", "X","M", "X","M", "X","M", "X", "M", "D", "M", "X", "M", "X", "M", "X", "M"]

        mIndex = 0
        cIndex = 0
        
        while mIndex < len(mOp) or cIndex < len(cOp):
            # If the current CIGAR operation is S, H, or I, add that to the output. The MD tag doesn't have these
            print mIndex
            print cIndex
            

            if cOp[cIndex] == "H" or cOp[cIndex] == "S" or cOp[cIndex] == "I":
                mergeOps.append(cOp[cIndex])
                mergeCts.append(cCt[cIndex])
                cIndex += 1

            # Special case: Insertions only occur in CIGAR. If an insertion follows a match, it is necessary to go back one step and shorten the match string.
            #if cOp[cIndex] == "I":
            #    print mergeCts[-1]
            #    exit()
            #    mergeCts[-1] = mergeCts[-1] - cCt[cIndex]
            #    mergeOps.append(cOp[cIndex])
            #    mergeCts.append(cCt[cIndex])
            #    cIndex += 1
 
 
            # Select the "shorter" operation
            else:
                if cCt[cIndex] < mCt[mIndex]:
                # If the CIGAR string lists fewer matched bases than MD, it means the CIGAR has an insertion not listed in MD
                    mCt[mIndex] = mCt[mIndex] - cCt[cIndex]
                    mergeOps.append(cOp[cIndex])
                    mergeCts.append(cCt[cIndex])
                    cIndex += 1

                if cCt[cIndex] > mCt[mIndex]:
                # If the CIGAR string lists more matched bases than MD, it means that MD has a mismatch not listed in CIGAR
                    cCt[cIndex] = cCt[cIndex] - mCt[mIndex]
                    mergeOps.append(mOp[mIndex])
                    mergeCts.append(mCt[mIndex])
                    mIndex += 1
                    

                if cCt[cIndex] == mCt[mIndex] and cOp[cIndex] == mOp[mIndex]:
                    mergeOps.append(mOp[mIndex])
                    mergeCts.append(mCt[mIndex])
                    mIndex += 1
                    cIndex += 1
            #if cCt[cIndex] == 0:
            #    cIndex += 1
            #if mCt[mIndex] == 0:
            #    mIndex += 1
            print mergeOps
            print mergeCts
            print "---------"
        exit()         


    def parseSpliceJunctions(self, genome):
        # This function takes the splice junction information from the SAM input and creates a SpliceJunction object for each.

        spliceJns = ((self.jM).split(":")[-1]).split(",")[1:]
        intronBounds = ((self.jI).split(":")[-1]).split(",")[1:]

        count = 0
        jnNum = 0
        jnObjects = [] 
        for entry in spliceJns:
            start = int(intronBounds[count])
            end = int(intronBounds[count + 1])
            sj = SpliceJunction(self.QNAME, jnNum, self.CHROM, start, end, self.strand, entry, genome)
            jnObjects.append(sj)

            # Check if junction is canonical or not. 
            if sj.isCanonical == False: self.isCanonical = False
            count += 2
            jnNum += 1
        
        return jnObjects

    def printableSAM(self, genome):
        # Returns a SAM-formatted string representation of the transcript
        if len(self.spliceJunctions) > 0:
            self.jI = "jI:B:i," + ",".join(str(i.pos) for i in self.getAllIntronBounds())
            self.jM = "jM:B:c," + ",".join(str(i) for i in self.getAllSJMotifs(genome))
        self.NM, self.MD = self.getNMandMDFlags(genome)        

        fields = [ self.QNAME, self.FLAG, self.CHROM, self.POS, self.MAPQ, self.CIGAR, self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, self.QUAL, self.otherFields, "NM:i:" + str(self.NM), self.MD, self.jM, self.jI ]
        return "\t".join([str(x) for x in fields]).strip()

    def printableFa(self):
        # Returns a fasta-formatted string representation of the transcript
        fasta1 = ">" + self.QNAME
        fastaSeq = [self.SEQ[i:i+80] for i in range(0, len(self.SEQ), 80)]
        return fasta1 + "\n" + "\n".join(fastaSeq)

    def getAllIntronBounds(self):
        # Return all intron bound objects belonging to this transcript

        result = []
        for jn in self.spliceJunctions:
            b = jn.bounds
            result.append(b[0])
            result.append(b[1])
        return result
   
    def getAllSJMotifs(self, genome):
    #    # Return all splice junction motifs translated into their numeric STAR codes
        result = []
        for jn in self.spliceJunctions:
            SpliceJunction.recheckJnStr(jn, genome)
            result.append(jn.jnStr)
        return result
 
    def getNMandMDFlags(self, genome):
        # This function uses the transcript sequence, its CIGAR string, and the reference genome to create NM and MD sam flags.
        NM = 0
        MD = "MD:Z:"
        MVal = 0
        seqPos = 0
        genomePos = self.POS

        operations, counts = self.splitCIGAR()
        for op, ct in zip(operations, counts):
            if op == "M":
                for i in range(0,ct):
                    currBase = self.SEQ[seqPos]
                    refBase = genome.sequence({'chr': self.CHROM, 'start': genomePos, 'stop': genomePos}, one_based=True).upper() 
                    if currBase != refBase:
                        # End any match we have going and add the mismatch
                        MD = MD + str(MVal) #if MVal > 0: 
                        MVal = 0 
                        MD = MD + refBase 
                        NM += 1
                    else:
                        MVal += 1
                    seqPos += 1
                    genomePos += 1
            if op == "D":
                # End any match we have going and add the missing reference bases
                MD = MD + str(MVal) #if MVal > 0: 
                MVal = 0
                refBases = genome.sequence({'chr': self.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True).upper() 
                MD = MD + "^" + refBases
                NM += ct
                genomePos += ct
            # For insertions and soft clips, we move on without adding to the MD
            if op in ["I", "S"]:
                seqPos += ct
                if op == "I": NM += ct
            if op in ["N", "H"]:
                genomePos += ct
                
        if MVal > 0: MD = MD + str(MVal) 
        return str(NM), MD
                  
