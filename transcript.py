# This file contains classes for the clean_splice_jns.py program

from spliceJunction import SpliceJunction

class Transcript:

    def __init__(self, sam):
        samFields = sam.strip().split('\t')

        # These attributes are initialized directly from the input SAM entry 
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
        self.NH = samFields[11]
        self.HI = samFields[12]
        self.NM = samFields[13]
        self.MD = samFields[14]
        self.jM = samFields[15]
        self.jI = samFields[16]
 
        # These attributes are set by parsing functions
        self.spliceJunctions = []
        self.isCanonical = True

        # Only run this section if there are splice junctions
        if "-1" not in self.jM:
            # Create an object for each splice junction
            self.spliceJunctions = self.parseSpliceJunctions()            


    def parseSpliceJunctions(self):
        # This function takes the splice junction information from the SAM input and creates a SpliceJunction object for each.

        spliceJns = ((self.jM).split(":")[-1]).split(",")[1:]
        intronBounds = ((self.jI).split(":")[-1]).split(",")[1:]

        count = 0
        jnObjects = [] 
        for entry in spliceJns:
            start = int(intronBounds[count])
            end = int(intronBounds[count + 1])
            sj = SpliceJunction(self.CHROM, start, end, entry)
            jnObjects.append(sj)

            # Check if junction is canonical or not. 
            if sj.isCanonical == False: self.isCanonical = False
            count += 2
        
        return jnObjects

    def printableSAM(self):
        # Returns a SAM-formatted string representation of the transcript
        fields = [ self.QNAME, self.FLAG, self.CHROM, self.POS, self.MAPQ, self.CIGAR, self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, self.QUAL, self.NH, self.HI, self.NM, self.MD, self.jM, self.jI ]
        return "\t".join([str(x) for x in fields])

    def printSpliceJunctions(self):
        
        for sj in self.spliceJunctions:
            print sj.getBED()

 
