# This file contains classes for the clean_splice_jns.py program

from spliceJunction import SpliceJunction
from intronBound import IntronBound
import pyfasta
import pybedtools
import re

class Transcript:

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
        #self.NH = ""
        #self.HI = ""
        self.NM = ""
        self.MD = ""
        self.jM = ""
        self.jI = ""        
        otherFields = []
        for field in samFields[11:len(samFields)]:
            #if field.startswith("NH"): self.NH = field
            #elif field.startswith("HI"): self.HI = field
            if field.startswith("NM"): self.NM = field
            elif field.startswith("MD"): self.MD = field
            elif field.startswith("jM"): self.jM = field 
            elif field.startswith("jI"): self.jI = field
            else: otherFields.append(field)

        self.otherFields = "\t".join(otherFields)        
        # These attributes are set by parsing the inputs
        self.spliceJunctions = []
        self.isCanonical = True
        self.strand = "+"        
        if self.FLAG == 16: self.strand = "-"

        #print self.SEQ
        #print self.referenceSeq

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

    #def getReferenceSequence(self, genome):
        # This function extracts the reference sequence of the region that the transcript mapped to. It uses POS, the start of 
        # the mapping, and gets the length of the mapping by summing the numbers in the CIGAR string

    #    alignTypes, counts = self.splitCIGAR()
    #    matchLen = sum(counts)
    #    seqStart = self.POS - 1 # Convert to 0-based
    #    seqEnd = seqStart + matchLen

    #    return genome[self.CHROM][seqStart:seqEnd].upper()

    def splitCIGAR(self):
        # Takes CIGAR string from SAM and splits it into two lists: one with capital letters (match operators), and one with the number of bases

        alignTypes = re.sub('[0-9]', " ", self.CIGAR).split()
        counts = re.sub('[A-Z]', " ", self.CIGAR).split()
        counts = [int(i) for i in counts]

        return alignTypes, counts


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
                    refBase = genome.sequence({'chr': self.CHROM, 'start': genomePos, 'stop': genomePos}, one_based=True)
                    if currBase.upper() != refBase.upper():
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
                refBases = genome.sequence({'chr': self.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True) 
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
                  
