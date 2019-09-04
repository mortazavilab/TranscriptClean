import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestDeletionCorr(object):

    def test_correctable_deletion(self):
        """ Toy transcript with sequence AA-GA, where the '-' is a deletion of 
            the base 'A'.
            chr1: 202,892,094 - 202,892,098. Deletion is at 202,892,096 """

        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "2M1D2M", "*",
                      "0", "0", "AAGA", "*",	"NM:i:1", "MD:Z:4", "jI:B:i,-1",
                      "jM:B:c,-1" ])

        genome = Fasta("input_files/hg38_chr1.fa")
        maxLen = 5
        spliceAnnot = None
        variants = {}
        transcriptErrorLog = open("scratch/TE.log", 'w')
        logInfo = TC.init_log_info()

        # Init transcript object
        transcript = t2.Transcript2(sam_fields, genome, spliceAnnot)

        # Run correction
        TC.correctDeletions(transcript, genome, variants, maxLen, logInfo,
                             transcriptErrorLog)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAAGA"
        assert transcript.CIGAR == "5M"

        transcriptErrorLog.close()

    def test_not_correctable_deletion(self):
        """ Same deletion again, but correction cutoff set to 0 """

        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "2M1D2M", "*",
                      "0", "0", "AAGA", "*",    "NM:i:1", "MD:Z:4", "jI:B:i,-1",
                      "jM:B:c,-1" ])

        genome = Fasta("input_files/hg38_chr1.fa")
        maxLen = 0
        spliceAnnot = None
        variants = {}
        transcriptErrorLog = open("scratch/TE.log", 'w')
        logInfo = TC.init_log_info()

        # Init transcript object
        transcript = t2.Transcript2(sam_fields, genome, spliceAnnot)

        # Run correction
        TC.correctDeletions(transcript, genome, variants, maxLen, logInfo,
                             transcriptErrorLog)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAGA"
        assert transcript.CIGAR == "2M1D2M"

        transcriptErrorLog.close()

    def test_variant_deletion(self):
        """ Same deletion again, but with a matching variant at the same 
            location. Correct action is to leave the deletion in place """

        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "2M1D2M", "*",
                      "0", "0", "AAGA", "*",    "NM:i:1", "MD:Z:4", "jI:B:i,-1",
                      "jM:B:c,-1" ])

        genome = Fasta("input_files/hg38_chr1.fa")
        maxLen = 5
        spliceAnnot = None
        variants = {"chr1_202892096_202892096": 1}
        transcriptErrorLog = open("scratch/TE.log", 'w')
        logInfo = TC.init_log_info()

        # Init transcript object
        transcript = t2.Transcript2(sam_fields, genome, spliceAnnot)

        # Run correction
        TC.correctDeletions(transcript, genome, variants, maxLen, logInfo,
                             transcriptErrorLog)

        # Check to see if deletion is still there as expected
        assert transcript.SEQ == "AAGA"
        assert transcript.CIGAR == "2M1D2M"

        transcriptErrorLog.close()

