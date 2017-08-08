# This file contains classes for the clean_splice_jns.py program

class SpliceJunction:

    def __init__(self, chrom, start, end, jnStr):
        
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        if int(jnStr) == 0 or jnStr == 20:
            self.isCanonical = False
        else:
            self.isCanonical = True   

    def isCanonical(self):
        return self.isCanonical

    def getBED(self):
        # Format the splice junction as a BED entry, with 0-based start and end. For now it is just a string, but later I might use an object
        return "\t".join([ self.chrom, str(self.start - 1), str(self.end) ]) 
