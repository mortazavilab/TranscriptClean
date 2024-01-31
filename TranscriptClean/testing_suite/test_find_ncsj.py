import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import transcript as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestFindNCSJ(object):

    def test_mark_canonical(self):
        """ In this test, we check whether TC correctly detects the junction
            and labels it canonical.

            Toy transcript with sequence AAG|GAA, where the splice motif (GT-AG)
            is canonical.
            chr1: 23,071,357 - 23,072,126 """
        test_dir = os.path.dirname(__file__)

        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "3M764N3M", "*",
                      "0", "0", "AAGGAA", "*",  "NM:i:0", "MD:Z:6"]

        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        maxLen = 5
        spliceAnnot = {}
        variants = {}

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Check if the intron bounds are correct
        intronBounds = transcript.getAllIntronBounds()
        assert intronBounds[0].pos == 23071360
        assert intronBounds[1].pos == 23072123

        # Check if the overall junction is labeled correctly
        assert (transcript.spliceJunctions[0]).isCanonical == True

        # Check if the overall transcript is labeled correctly
        assert transcript.isCanonical == True

    def test_mark_noncanonical(self):
        """ In this test, we check whether TC correctly detects the junction
            and labels it noncanonical.

            Toy transcript with sequence GGT|GTG, where the splice motif (AA-CA)
            is noncanonical.
            chr1: 23,072,197 - 23,073,291.  """
        test_dir = os.path.dirname(__file__)

        sam_fields = ["test_read", "0", "chr1", "23072197", "255", "3M1091N3M", "*",
                      "0", "0", "GGTGTG", "*",	"NM:i:0", "MD:Z:6"]

        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        maxLen = 5
        spliceAnnot = {}
        variants = {}

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Check if the intron bounds are correct
        intronBounds = transcript.getAllIntronBounds()
        assert intronBounds[0].pos == 23072200
        assert intronBounds[1].pos == 23073290

        # Check if the overall junction is labeled correctly
        assert (transcript.spliceJunctions[0]).isCanonical == False

        # Check if the overall transcript is labeled correctly
        assert transcript.isCanonical == False
