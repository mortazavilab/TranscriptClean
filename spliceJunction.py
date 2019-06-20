# This file contains the SpliceJunction class and associated functions

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

    #def isCanonical(self):
    #    """ If both intron bounds of the junction are canonical, then so is the 
    #        splice junction as a whole."""
    #            
    #    return self.bounds[0].isCanonical and self.bounds[1].isCanonical

    def get_splice_donor(self):
        """ Return the IntronBound object that is the splice donor """
        if self.strand == "+":
            return self.bounds[0]
        elif self.strand == "-":
            return self.bounds[1]

    def get_splice_acceptor(self):
        """ Return the IntronBound object that is the splice acceptor """
        if self.strand == "+":
            return self.bounds[1]
        elif self.strand == "-":
            return self.bounds[0]

    def recheckPosition(self):
        """ Get start and end position from its intron bounds """
        self.start = self.bounds[0].pos
        self.end = self.bounds[1].pos

    def recheckJnStr(self, genome, spliceAnnot):
        """ Check the splice junction sequence motif to determine whether the 
            jnStr should be changed"""

        startBases = IntronBound.getSpliceMotif(self.bounds[0], genome) 
        endBases = IntronBound.getSpliceMotif(self.bounds[1], genome).upper()        

        if (self.chrom + "_" + str(self.start)) in spliceAnnot and (self.chrom + "_" + str(self.end)) in spliceAnnot:
            motifCode = 20 
        else: motifCode = 0

        motifCode += getSJMotifCode(startBases, endBases)
        self.jnStr = str(motifCode)
        self.isCanonical = self.jnStr != "0" and self.jnStr != "20" #self.bounds[0].isCanonical and self.bounds[1].isCanonical
        return        
        
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
