import pytest
from pyfaidx import Fasta
import sys
sys.path.append("..")
import transcript as ts
import TranscriptClean as TC
@pytest.mark.unit

class TestPrintFaSeq(object):
    def test_plus_strand(self):
        """ Toy transcript with sequence AAAGA on the plus strand. Sequence
            should be output as listed in the SAM fields."""

        sam_fields = ["test_read", "0", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AAAGA", "*",	"NM:i:0", "MD:Z:5", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None

        # Init transcript object
        transcript = ts.Transcript(sam_fields, genome, spliceAnnot)

        # Output fasta and check against expected
        expected_fasta = ">test_read" + "\n" + "AAAGA"

        assert transcript.printableFa() == expected_fasta

    def test_minus_strand(self):
        """ Toy transcript on the minus strand with SAM sequence AAAGA.
            Sequence for FASTA file should be the reverse complement.
        """

        sam_fields = ["test_read", "16", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AAAGA", "*",   "NM:i:0", "MD:Z:5", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None

        # Init transcript object
        transcript = ts.Transcript(sam_fields, genome, spliceAnnot)

        # Output fasta and check against expected
        expected_fasta = ">test_read" + "\n" + "TCTTT"

        assert transcript.printableFa() == expected_fasta
