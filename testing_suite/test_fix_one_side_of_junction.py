import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import spliceJunction as sj
import intronBound as ib
import TranscriptClean as TC

class TestFixSideOfJunction(object):

    def test_fix_donor_plus(self):
        """ Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126
        """

        # Process references
        sjFile = "input_files/test_junctions.txt"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)
        genome = Fasta("input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = "\t".join(["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"])
        transcript = t2.Transcript2(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM, 
                                                         transcript.POS, jnNumber, 
                                                         donor, 2, genome, 
                                                         transcript.SEQ, 
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"


    def test_fix_acceptor_plus(self):
        """ Toy transcript with sequence AAG|AA, where the splice motif
            is noncanonical but located 1 bp from a canonical splice acceptor.
            chr1: 23,071,357 - 23,072,126
        """ 
        # Process references
        sjFile = "input_files/test_junctions.txt"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)
        genome = Fasta("input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = "\t".join(["test_read", "0", "chr1", "23071357", "255", "3M765N2M", "*",
                      "0", "0", "AAGAA", "*",  "NM:i:0", "MD:Z:5"])
        transcript = t2.Transcript2(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        acceptor = (transcript.spliceJunctions[jnNumber]).bounds[1]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM,
                                                         transcript.POS, jnNumber,
                                                         acceptor, -1, genome,
                                                         transcript.SEQ,
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"
