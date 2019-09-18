import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.unit

class TestInitTranscript(object):
    def test_unmapped_read(self):
        """ The supplied read is unmapped. This means that a limited object is 
            created, and the logInfo struct logs the unmapped status."""

        sam_file = "input_files/init_transcript/unmapped.sam"
        with open(sam_file, 'r') as f:
            sam_line = f.readline().strip()

        genome = None
        sjAnnot = None
      
        transcript, logInfo = TC.transcript_init(sam_line, genome, sjAnnot)
        assert transcript == None
        assert logInfo.mapping == "unmapped"
        assert logInfo.corrected_deletions == \
               logInfo.uncorrected_deletions == \
               logInfo.variant_deletions == \
               logInfo.corrected_insertions == \
               logInfo.uncorrected_insertions == \
               logInfo.variant_insertions == \
               logInfo.corrected_mismatches == \
               logInfo.uncorrected_mismatches == \
               logInfo.corrected_NC_SJs == logInfo.uncorrected_NC_SJs == "NA"
