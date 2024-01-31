import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import transcript as t2
import TranscriptClean.spliceJunction as sj
import TranscriptClean.intronBound as ib
import TranscriptClean.TranscriptClean as TC

class TestFixSideOfJunction(object):

    def test_fix_donor_case1(self):
        """ Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

            So-called case # 1
        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"]
        transcript = t2.Transcript(sam_fields, genome, sjDict)
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


    def test_fix_acceptor_case2(self):
        """ Toy transcript with sequence AAG|AA, where the splice motif
            is noncanonical but located 1 bp from a canonical splice acceptor.
            chr1: 23,071,357 - 23,072,126

            So-called case #2
        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "3M765N2M", "*",
                      "0", "0", "AAGAA", "*",  "NM:i:0", "MD:Z:5"]
        transcript = t2.Transcript(sam_fields, genome, sjDict)
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

    def test_fix_donor_case3(self):
        """ Toy transcript with sequence AAGGT|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

            So-called case #3
        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "5M762N3M", "*",
                      "0", "0", "AAGGTGAA", "*",  "NM:i:0", "MD:Z:8"]
        transcript = t2.Transcript(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM,
                                                         transcript.POS, jnNumber,
                                                         donor, -2, genome,
                                                         transcript.SEQ,
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"

    def test_fix_acceptor_case4(self):
        """ Toy transcript with sequence AAG|GGAA, where the splice motif
            is noncanonical but located 1 bp from a canonical splice acceptor.
            chr1: 23,071,357 - 23,072,126

            So-called case #4
        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "3M763N4M", "*",
                      "0", "0", "AAGGGAA", "*",  "NM:i:0", "MD:Z:7"]
        transcript = t2.Transcript(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        acceptor = (transcript.spliceJunctions[jnNumber]).bounds[1]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM,
                                                         transcript.POS, jnNumber,
                                                         acceptor, 1, genome,
                                                         transcript.SEQ,
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"
