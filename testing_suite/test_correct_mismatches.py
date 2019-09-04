import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestMismatchCorr(object):

    def test_correctable_mismatch(self):
        """ Toy transcript with sequence AACGA, where the C is a mismatch to the
            reference base 'A'.
            chr1: 202,892,094 - 202,892,098. Mismatch is at 202,892,096 """

        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AACGA", "*",	"NM:i:1", "MD:Z:2A2", "jI:B:i,-1",
                      "jM:B:c,-1" ])

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None
        variants = {}
        transcriptErrorLog = open("scratch/TE.log", 'w')
        logInfo = TC.init_log_info()

        # Init transcript object
        transcript = t2.Transcript2(sam_fields, genome, spliceAnnot)

        # Run correction
        TC.correctMismatches(transcript, genome, variants, logInfo,
                             transcriptErrorLog)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAAGA"
        assert transcript.CIGAR == "5M"

        transcriptErrorLog.close()

    def test_variant_mismatch(self):
        """ Toy transcript with sequence AACGA, where the C is a mismatch to the
            reference base 'A', but is a known SNP.
            chr1: 202,892,094 - 202,892,098. Mismatch is at 202,892,096 """
        
        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AACGA", "*",   "NM:i:1", "MD:Z:2A2", "jI:B:i,-1",
                      "jM:B:c,-1" ])

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None
        variants = {"chr1_202892096" : ["C", "T"] }
        transcriptErrorLog = open("scratch/TE.log", 'w')
        logInfo = TC.init_log_info()

        # Init transcript object
        transcript = t2.Transcript2(sam_fields, genome, spliceAnnot)

        # Run correction
        TC.correctMismatches(transcript, genome, variants, logInfo,
                             transcriptErrorLog)

        # Check to see if correction was successful
        assert transcript.SEQ == "AACGA"
        assert transcript.CIGAR == "5M"

        transcriptErrorLog.close()

    def test_wrong_variant_mismatch(self):
        """ Toy transcript with sequence AACGA, where the C is a mismatch to the
            reference base 'A' in the location, but not matching, a known SNP.
            chr1: 202,892,094 - 202,892,098. Mismatch is at 202,892,096 """

        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "5M", "*",
                      "0", "0", "AACGA", "*",   "NM:i:1", "MD:Z:2A2", "jI:B:i,-1",
                      "jM:B:c,-1" ])

        genome = Fasta("input_files/hg38_chr1.fa")
        spliceAnnot = None
        variants = {"chr1_202892096" : ["G"] }
        transcriptErrorLog = open("scratch/TE.log", 'w')
        logInfo = TC.init_log_info()

        # Init transcript object
        transcript = t2.Transcript2(sam_fields, genome, spliceAnnot)

        # Run correction
        TC.correctMismatches(transcript, genome, variants, logInfo,
                             transcriptErrorLog)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAAGA"
        assert transcript.CIGAR == "5M"

        transcriptErrorLog.close()
