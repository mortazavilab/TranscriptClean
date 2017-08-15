# This file contains classes for the clean_splice_jns.py program

from spliceJunction import SpliceJunction
from intronBound import IntronBound
import pyfasta
import pybedtools
import re

class Transcript:

    def __init__(self, sam, genome):
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
        self.jM = samFields[-2]
        self.jI = samFields[-1]
 
        # These attributes are set by parsing the inputs
        self.spliceJunctions = []
        self.isCanonical = True
        self.strand = "+"
        self.referenceSeq = self.getReferenceSequence(genome)
        if self.FLAG == 16: self.strand = "-"

        #print self.SEQ
        #print self.referenceSeq

        # Only run this section if there are splice junctions
        if "-1" not in self.jM:
            # Create an object for each splice junction
            self.spliceJunctions = self.parseSpliceJunctions()            
 
    def recheckCanonical(self):
        for jn in self.spliceJunctions:
            if jn.isCanonical == False:
                self.isCanonical = False
                return False
        self.isCanonical = True
        return True

    def getReferenceSequence(self, genome):
        # This function extracts the reference sequence of the region that the transcript mapped to. It uses POS, the start of 
        # the mapping, and gets the length of the mapping by summing the numbers in the CIGAR string

        alignTypes, counts = self.splitCIGAR()
        matchLen = sum(counts)
        seqStart = self.POS - 1 # Convert to 0-based
        seqEnd = seqStart + matchLen

        return genome[self.CHROM][seqStart:seqEnd].upper()

    def splitCIGAR(self):
        # Takes CIGAR string from SAM and splits it into two lists: one with capital letters (match operators), and one with the number of bases
        # The relative order of the elements is maintained, but the lists are reversed for cases on the '-' strand    

        alignTypes = re.sub('[0-9]', " ", self.CIGAR).split()
        counts = re.sub('[A-Z]', " ", self.CIGAR).split()
        counts = [int(i) for i in counts]

        return alignTypes, counts


    def parseSpliceJunctions(self):
        # This function takes the splice junction information from the SAM input and creates a SpliceJunction object for each.

        spliceJns = ((self.jM).split(":")[-1]).split(",")[1:]
        intronBounds = ((self.jI).split(":")[-1]).split(",")[1:]

        count = 0
        jnObjects = [] 
        for entry in spliceJns:
            start = int(intronBounds[count])
            end = int(intronBounds[count + 1])
            sj = SpliceJunction(self.QNAME, count, self.CHROM, start, end, self.strand, entry)
            jnObjects.append(sj)

            # Check if junction is canonical or not. 
            if sj.isCanonical == False: self.isCanonical = False
            count += 2
        
        return jnObjects

    def printableSAM(self):
        # Returns a SAM-formatted string representation of the transcript
        fields = [ self.QNAME, self.FLAG, self.CHROM, self.POS, self.MAPQ, self.CIGAR, self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, self.QUAL, self.jM, self.jI ]
        return "\t".join([str(x) for x in fields])


    def getAllIntronBounds(self):
        # Return all intron bound objects belonging to this transcript

        result = []
        for jn in self.spliceJunctions:
            b = jn.bounds
            result.append(b[0])
            result.append(b[1])
        return result
    
 
