import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import transcriptclean.transcript as t2
import transcriptclean.TranscriptClean as TC
@pytest.mark.unit

class TestComputeMD(object):
    def test_perfect_match_no_introns(self):
        """ Compute the correct MD tag for a transcript that is a perfect
            reference match with no introns. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        correct_MD = "MD:Z:155"
        correct_NM = "NM:i:0"
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM

    def test_perfect_match_with_introns(self):
        """ Compute the correct MD tag for a transcript that is a perfect
            reference match containing introns. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_twoIntrons.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        correct_MD = "MD:Z:3400"
        correct_NM = "NM:i:0"
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM

    def test_mismatch(self):
        """ Compute the correct MD tag for a spliced transcript that contains a
            mismatch. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/mismatch.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        correct_MD = "MD:Z:2G2813"
        correct_NM = "NM:i:1"
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM

    def test_deletion(self):
        """ Compute the correct MD tag for a spliced transcript that contains a
            deletion. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/deletion.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        correct_MD = "MD:Z:2629^G377"
        correct_NM = "NM:i:1"
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM

    def test_insertion(self):
        """ Compute the correct MD tag for a spliced transcript that contains an
            insertion. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/insertion.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        correct_MD = "MD:Z:3069"
        correct_NM = "NM:i:2"
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM

    def test_deletion_insertion_mismatch(self):
        """ Compute the correct MD tag for a spliced transcript that contains an
            insertion, deletion, and mismatch. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/deletion_insertion_mismatch.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        correct_MD = "MD:Z:475C0A0C0C0A1082^C1347G17C205"
        correct_NM = "NM:i:9"
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM

    def test_insertion_deletion_mismatch_ncsj(self):
        """ Compute the correct MD tag for a transcript that contains an
           insertion, deletion, mismatch, and noncanonical splice junction in
           it. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/deletion_insertion_mismatch_nc.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        correct_MD = "MD:Z:414G0A450^C405^C2435"
        correct_NM = "NM:i:5"
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM
