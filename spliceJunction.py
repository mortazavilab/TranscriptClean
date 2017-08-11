# This file contains classes for the clean_splice_jns.py program

import pybedtools

class SpliceJunction:

    def __init__(self, chrom, start, end, strand, jnStr):
        
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        if int(jnStr) == 0:
            self.isCanonical = False
        else:
            self.isCanonical = True   

    def isCanonical(self):
        return self.isCanonical

    def getBED(self, mode):
        # Format the splice junction as a BedTool object, with 0-based start and end. 
        # If mode == "start", we return a bed for the start position.
        # If mode == "end", we return a bed for the end position.
        if mode == "start": pos = self.start
        if mode == "end": pos = self.end
        bedstr = "\t".join([ self.chrom, str(pos - 1), str(pos), ".", "0", self.strand ])
        return pybedtools.BedTool(bedstr, from_string= True) 
