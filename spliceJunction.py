# This file contains classes for the clean_splice_jns.py program

from intronBound import IntronBound

class SpliceJunction:

    def __init__(self, transcriptID, jnNumber, chrom, start, end, strand, jnStr):
        
        self.transcriptID = transcriptID
        self.jnNumber = int(jnNumber)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        if int(jnStr) == 0:
            self.isCanonical = False
        else:
            self.isCanonical = True   

        # Create an intronBound object for each end of the junction
        left = IntronBound(self.transcriptID, self.jnNumber, "0", self.chrom, self.start, self.strand, jnStr)
        right = IntronBound(self.transcriptID, self.jnNumber, "1", self.chrom, self.end, self.strand, jnStr)
        self.bounds = [left, right]

    def isCanonical(self):
        # If both intron bounds of the junction are canonical, then so is the splice juntion as a whole.
                
        return self.bounds[0].isCanonical and self.bounds[1].isCanonical


