import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestMismatchCorr(object):

    def test_correctable_mismatch(self):
        """ Toy transcript with sequence AACGA, where the C is a mismatch to the
            reference base 'A'.
            chr1: 202,892,094 - 202,892,098. Mismatch is at 202,892,096 """

        sam_fields = ["test_read", "0", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AACGA", "*",	"NM:i:1", "MD:Z:2A2", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None
        variants = {}
        logInfo = TC.init_log_info(sam_fields)

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Run correction
        error_entries = TC.correctMismatches(transcript, genome, variants, logInfo)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAAGA"
        assert transcript.CIGAR == "5M"

        # Check the number and content of the transcript error entries
        assert len(error_entries) == 1
        assert "Corrected" in error_entries[0]


    def test_variant_mismatch(self):
        """ Toy transcript with sequence AACGA, where the C is a mismatch to the
            reference base 'A', but is a known SNP.
            chr1: 202,892,094 - 202,892,098. Mismatch is at 202,892,096 """
        
        sam_fields = ["test_read", "0", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AACGA", "*",   "NM:i:1", "MD:Z:2A2", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None
        variants = {"chr1_202892096" : ["C", "T"] }
        logInfo = TC.init_log_info(sam_fields)

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Run correction
        error_entries = TC.correctMismatches(transcript, genome, variants, logInfo)

        # Check to see if correction was successful
        assert transcript.SEQ == "AACGA"
        assert transcript.CIGAR == "5M"

        # Check the number and content of the transcript error entries
        assert len(error_entries) == 1
        assert "Uncorrected" in error_entries[0]
        assert "VariantMatch" in error_entries[0]

    def test_wrong_variant_mismatch(self):
        """ Toy transcript with sequence AACGA, where the C is a mismatch to the
            reference base 'A' in the location, but not matching, a known SNP.
            chr1: 202,892,094 - 202,892,098. Mismatch is at 202,892,096 """

        sam_fields = ["test_read", "0", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AACGA", "*",   "NM:i:1", "MD:Z:2A2", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None
        variants = {"chr1_202892096" : ["G"] }
        logInfo = TC.init_log_info(sam_fields)

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Run correction
        error_entries = TC.correctMismatches(transcript, genome, variants, logInfo)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAAGA"
        assert transcript.CIGAR == "5M"

        # Check the number and content of the transcript error entries
        assert len(error_entries) == 1
        assert "Corrected" in error_entries[0]
        assert "VariantMatch" not in error_entries[0]

