import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import transcriptclean.transcript as t2
import transcriptclean.TranscriptClean as TC
@pytest.mark.unit

class TestInsertionCorr(object):

    def test_correctable_insertion(self):
        """ Toy transcript with sequence AAATTGA, where the Ts are a 2 bp insertion.
            chr1: 202,892,094 - 202,892,098. Insertion is between position
            202,892,096 and 202,892,097. The genomic position used to refer
            to it is 202,892,097"""
        test_dir = os.path.dirname(__file__)

        sam_fields = ["test_read", "0", "chr1", "202892094", "255", "3M2I2M", "*",
                      "0", "0", "AAATTGA", "*",	"NM:i:2", "MD:Z:5", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        maxLen = 5
        spliceAnnot = None
        variants = {}
        logInfo = TC.init_log_info(sam_fields)

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Run correction
        TE_entries = TC.correctInsertions(transcript, genome, variants, maxLen, logInfo)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAAGA"
        assert transcript.CIGAR == "5M"

        # Check the log entries
        expected_log = "\t".join(["test_read", "chr1_202892096_202892098",
                                  "Insertion", "2", "Corrected", "NA"]) + "\n"
        assert TE_entries == expected_log


    def test_not_correctable_insertion(self):
        """ Toy transcript with sequence AAATTGA, where the Ts are a 2 bp insertion.
            chr1: 202,892,094 - 202,892,098. Insertion is between position
            202,892,096 and 202,892,097. The genomic position used to refer
            to it is 202,892,097"""
        test_dir = os.path.dirname(__file__)

        sam_fields = ["test_read", "0", "chr1", "202892094", "255", "3M2I2M", "*",
                      "0", "0", "AAATTGA", "*", "NM:i:2", "MD:Z:5", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        maxLen = 1
        spliceAnnot = None
        variants = {}
        logInfo = TC.init_log_info(sam_fields)

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Run correction
        TE_entries = TC.correctInsertions(transcript, genome, variants, maxLen, logInfo)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAATTGA"
        assert transcript.CIGAR == "3M2I2M"

        # Check the log entries
        expected_log = "\t".join(["test_read", "chr1_202892096_202892098",
                                  "Insertion", "2", "Uncorrected", "TooLarge"]) + "\n"
        assert TE_entries == expected_log


    def test_variant_insertion(self):
        """ Toy transcript with sequence AAATTGA, where the Ts are a 2 bp
            insertion that matches a known variant.
            chr1: 202,892,094 - 202,892,098. Insertion is between position
            202,892,096 and 202,892,097. The genomic position used to refer
            to it is 202,892,097 """
        test_dir = os.path.dirname(__file__)

        sam_fields = ["test_read", "0", "chr1", "202892094", "255", "3M2I2M", "*",
                      "0", "0", "AAATTGA", "*", "NM:i:2", "MD:Z:5", "jI:B:i,-1",
                      "jM:B:c,-1" ]

        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        maxLen = 5
        spliceAnnot = None
        variants = {"chr1_202892096_202892098": "TT"}
        logInfo = TC.init_log_info(sam_fields)

        # Init transcript object
        transcript = t2.Transcript(sam_fields, genome, spliceAnnot)

        # Run correction
        TE_entries = TC.correctInsertions(transcript, genome, variants, maxLen, logInfo)

        # Check to see if correction was successful
        assert transcript.SEQ == "AAATTGA"
        assert transcript.CIGAR == "3M2I2M"

        # Check the log entries
        expected_log = "\t".join(["test_read", "chr1_202892096_202892098",
                                  "Insertion", "2", "Uncorrected", "VariantMatch"]) + "\n"
        assert TE_entries == expected_log
