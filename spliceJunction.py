# This file contains classes for the clean_splice_jns.py program

import pybedtools

class SpliceJunction:

    def __init__(self, transcriptID, jnNumber, chrom, start, end, strand, jnStr):
        
        self.transcriptID = transcriptID
        self.jnNumber = jnNumber
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

    def changeToCanonical(self):
        # This function is used when there is a canonical splice junction within 5 bp of this noncanonical junction.
        # It changes the start and end of this splice junction to those of the canonical one.
        if self.isCanonical == True: return
        
        # Convert new start and end to 1-based since they will be coming in as bed
        #self.start = newStart + 1
        #self.end = newEnd
        self.isCanonical = True
        return


    def getBED(self, mode):
        # Format the splice junction as a BedTool object, with 0-based start and end. 
        # If mode == "start", we return a bed for the start position.
        # If mode == "end", we return a bed for the end position.
        if mode == "start": 
            pos = self.start
            side = "0"
        if mode == "end": 
            pos = self.end
            side = "1"
        name = self.transcriptID + "__" + str(self.jnNumber) + "__" + side
        bedstr = "\t".join([ self.chrom, str(pos - 1), str(pos), name, "0", self.strand ])
        return bedstr #pybedtools.BedTool(bedstr, from_string= True) 
