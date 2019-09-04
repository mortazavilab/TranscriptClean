import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestComputeMD(object):
    def test_insertion_deletion_mismatch_ncsj(self):
        """ Compute the correct MD tag for a transcript that contains an 
           insertion, deletion, mismatch, and noncanonical splice junction in 
           it. """
 
        sam = "input_files/sams/deletion_insertion_mismatch_nc.sam"
        genome = Fasta("input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, {})
            transcript.getNMandMDFlags(genome)

        correct_MD = sam_line.split("\t")[14]
        assert transcript.MD == correct_MD
        
