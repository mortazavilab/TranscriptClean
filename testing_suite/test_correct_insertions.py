import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestInsertionCorr(object):

    def test_correctable_insertion(self):
        """ Toy transcript with sequence AAATTGA, where the Ts are a 2 bp insertion.
            chr1: 202,892,094 - 202,892,098. Insertion is at 202,892,096 """

        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "3M2I2M", "*",
                      "0", "0", "AAATTGA", "*",	"NM:i:2", "MD:Z:5", "jI:B:i,-1",
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
        TC.correctInsertions(transcript, genome, variants, maxLen, logInfo,
                             transcriptErrorLog)
    
        # Check to see if correction was successful
        assert transcript.SEQ == "AAAGA"
        assert transcript.CIGAR == "5M"

        transcriptErrorLog.close()

    def test_not_correctable_insertion(self):
        """ Toy transcript with sequence AAATTGA, where the Ts are a 2 bp insertion.
            chr1: 202,892,094 - 202,892,098. Insertion is at 202,892,096 """

        sam_fields = "\t".join(["test_read", "0", "chr1", "202892094", "255", "3M2I2M", "*",
                      "0", "0", "AAATTGA", "*", "NM:i:2", "MD:Z:5", "jI:B:i,-1",
                      "jM:B:c,-1" ])

        genome = Fasta("input_files/hg38_chr1.fa")
        maxLen = 1
        spliceAnnot = None
        variants = {}
        transcriptErrorLog = open("scratch/TE.log", 'w')
        logInfo = TC.init_log_info()

        # Init transcript object
        transcript = t2.Transcript2(sam_fields, genome, spliceAnnot)

        # Run correction
        TC.correctInsertions(transcript, genome, variants, maxLen, logInfo,
                             transcriptErrorLog)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAATTGA"
        assert transcript.CIGAR == "3M2I2M"

        transcriptErrorLog.close()
         
        
         
