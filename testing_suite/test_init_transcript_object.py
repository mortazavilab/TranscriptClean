import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.unit

class TestInitTranscript(object):
    def test_unmapped_read(self):
        """ The supplied read is unmapped. This means that no transcript object is 
            created, and the logInfo struct notes the unmapped status."""

        sam_file = "input_files/init_transcript/unmapped.sam"
        with open(sam_file, 'r') as f:
            sam_line = f.readline().strip()

        genome = None
        sjAnnot = set()
      
        transcript, logInfo = TC.transcript_init(sam_line, genome, sjAnnot)
        assert transcript == None
        assert logInfo.Mapping == "unmapped"
        assert logInfo.corrected_deletions == \
               logInfo.uncorrected_deletions == \
               logInfo.variant_deletions == \
               logInfo.corrected_insertions == \
               logInfo.uncorrected_insertions == \
               logInfo.variant_insertions == \
               logInfo.corrected_mismatches == \
               logInfo.uncorrected_mismatches == \
               logInfo.corrected_NC_SJs == logInfo.uncorrected_NC_SJs == "NA"

    def test_nonprimary_read(self):
        """ The supplied read is a non-primary alignment. This means that no 
            transcript object is created, and the logInfo struct notes the 
            non-primary status."""

        sam_file = "input_files/init_transcript/non-primary.sam"
        with open(sam_file, 'r') as f:
            sam_line = f.readline().strip()

        genome = None
        sjAnnot = set()

        transcript, logInfo = TC.transcript_init(sam_line, genome, sjAnnot)
        assert transcript == None
        assert logInfo.Mapping == "non-primary"
        assert logInfo.corrected_deletions == \
               logInfo.uncorrected_deletions == \
               logInfo.variant_deletions == \
               logInfo.corrected_insertions == \
               logInfo.uncorrected_insertions == \
               logInfo.variant_insertions == \
               logInfo.corrected_mismatches == \
               logInfo.uncorrected_mismatches == \
               logInfo.corrected_NC_SJs == logInfo.uncorrected_NC_SJs == "NA"

    def test_primary_monoexon_read(self):
        """ The supplied read is a primary alignment. This means that a
            transcript object is created, and the logInfo struct notes the
            primary status."""

        sam_file = "input_files/sams/perfectReferenceMatch_noIntrons.sam"
        with open(sam_file, 'r') as f:
            sam_line = f.readline().strip()

        genome = Fasta("input_files/hg38_chr1.fa")
        sjAnnot = set()

        transcript, logInfo = TC.transcript_init(sam_line, genome, sjAnnot)
        assert transcript.QNAME == "c21031/f2p3/3400"
        assert transcript.FLAG == 0
        assert transcript.CHROM == "chr1"
        assert transcript.POS == 192575775
        assert transcript.CIGAR == "155M"
        assert transcript.MD == "MD:Z:155"
        assert logInfo.Mapping == "primary"      
        assert logInfo.corrected_deletions == \
               logInfo.uncorrected_deletions == \
               logInfo.variant_deletions == \
               logInfo.corrected_insertions == \
               logInfo.uncorrected_insertions == \
               logInfo.variant_insertions == \
               logInfo.corrected_mismatches == \
               logInfo.uncorrected_mismatches == \
               logInfo.corrected_NC_SJs == logInfo.uncorrected_NC_SJs == "NA" 


