# This file contains the SpliceJunction class and associated functions

from intronBound import IntronBound

class SpliceJunction:

    def __init__(self, transcriptID, jnNumber, chrom, start, end, strand, genome, 
                 spliceAnnot):
        
        self.transcriptID = transcriptID
        self.jnNumber = int(jnNumber)
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

        # Create an intronBound object for each end of the junction
        left = IntronBound(self.transcriptID, self.jnNumber, "0", self.chrom,
                           self.start, self.strand)
        right = IntronBound(self.transcriptID, self.jnNumber, "1", self.chrom,
                            self.end, self.strand)
        self.bounds = [left, right]

        # Get splice motif
        #self.motif_code = self.isCanonical = None
        self.checkSpliceMotif(genome, spliceAnnot)
        
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

    def checkSpliceMotif(self, genome, spliceAnnot):
        """ Check the splice junction sequence motif to determine whether the 
            jnStr should be changed"""

        startBases = IntronBound.getSpliceMotif(self.bounds[0], genome) 
        endBases = IntronBound.getSpliceMotif(self.bounds[1], genome).upper()        

        # ID the splice donor/acceptor identity
        if self.strand == "+":
            type1 = "donor"
            type2 = "acceptor"
        elif self.strand == "-":
            type1 = "acceptor"
            type2 = "donor"
        if ("_".join([self.chrom, str(self.start), self.strand, type1])) in spliceAnnot and \
           ("_".join([self.chrom, str(self.end), self.strand, type2])) in spliceAnnot:
            motifCode = 20 
        else: 
            motifCode = 0

        motifCode += getSJMotifCode(startBases, endBases)
        self.motif_code = str(motifCode)
        self.isCanonical = self.motif_code != "0" and self.motif_code != "20"
                
        
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
