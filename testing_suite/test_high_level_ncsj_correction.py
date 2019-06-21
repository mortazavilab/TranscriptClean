import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import spliceJunction as sj
import intronBound as ib
import TranscriptClean as TC
import dstruct as dstruct

class TestNCSJCorrection(object):

   def test_correct_ncsj(self):
        """ Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

        """

        # Process references
        sjFile = "input_files/test_junctions.txt"
        outprefix = "scratch/test"
        TElog = open("scratch/test_clean.TE.log", 'w')
        refs = dstruct.Struct()
        refs.donors, refs.acceptors, refs.sjDict = TC.processSpliceAnnotation(sjFile, outprefix)
        refs.genome = Fasta("input_files/hg38_chr1.fa")

        # Init transcript object
        sam_fields = "\t".join(["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"])
        transcript = t2.Transcript2(sam_fields, refs.genome, refs.sjDict)
        jnNumber = 0
        maxDist = 5
        logInfo = TC.init_log_info()
        logInfo.TranscriptID = transcript.QNAME
        logInfo.Mapping = "primary"
        donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        assert transcript.isCanonical == False

        # Attempt to correct the splice junction
        transcript = TC.cleanNoncanonical(transcript, refs, maxDist, logInfo, TElog)

        assert transcript.isCanonical == True
        assert transcript.spliceJunctions[jnNumber].isCanonical == True
        assert transcript.SEQ == "AAGGAA"
        assert transcript.CIGAR == "3M764N3M"
        assert transcript.MD == "MD:Z:6"
        assert logInfo.corrected_NC_SJs == 1 
