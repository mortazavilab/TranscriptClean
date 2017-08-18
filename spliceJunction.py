# This file contains classes for the clean_splice_jns.py program

from intronBound import IntronBound

class SpliceJunction:

    def __init__(self, transcriptID, jnNumber, chrom, start, end, strand, jnStr, genome):
        
        self.transcriptID = transcriptID
        self.jnNumber = int(jnNumber)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.jnStr = jnStr
        if int(jnStr) == 0:
            self.isCanonical = False
        else:
            self.isCanonical = True   
        
        # Create an intronBound object for each end of the junction
        left = IntronBound(self.transcriptID, self.jnNumber, "0", self.chrom, self.start, self.strand, jnStr, genome)
        right = IntronBound(self.transcriptID, self.jnNumber, "1", self.chrom, self.end, self.strand, jnStr, genome)
        self.bounds = [left, right]

    def isCanonical(self):
        # If both intron bounds of the junction are canonical, then so is the splice juntion as a whole.
                
        return self.bounds[0].isCanonical and self.bounds[1].isCanonical

    def recheckJnStr(self, genome):
        # Check the splice junction sequence motif to determine whether the jnStr shouls be changed

        motif = IntronBound.getSpliceMotif(self.bounds[0], genome) + IntronBound.getSpliceMotif(self.bounds[1], genome)

        if motif == "GTAG":
            self.jnStr = "21"
        elif motif == "CTAC":
            self.jnStr = "22"
        elif motif == "GCAG":
            self.jnStr = "23"
        elif motif == "CTGC":
            self.jnStr = "24"
        elif motif == "ATAC":
            self.jnStr = "25"
        elif motif == "GTAT":
            self.jnStr = "26"
        else:
            self.jnStr = "0"
        return        
        

