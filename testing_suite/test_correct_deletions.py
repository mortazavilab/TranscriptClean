import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestDeletionCorr(object):

    def test_correctable_deletion(self):
        """ Toy transcript with sequence AAATTGA, where the Ts are a 2 bp insertion.
            Toy transcript with sequence AA-GA, where the '-' is a deletion of 
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
