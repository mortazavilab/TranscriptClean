import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import transcriptclean.transcript as t2
import transcriptclean.TranscriptClean as TC
@pytest.mark.unit

class TestSeqCigar(object):
    def test_perfect_match(self):
        """ Since this read is a perfect match to the reference, its CIGAR and
            sequence fields should definitely be the same length """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        sjDict = set()

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, sjDict)

        assert t2.check_seq_and_cigar_length(transcript.SEQ, transcript.CIGAR) == True

    def test_pre_correction_dmel(self):
        """ This is a noisy Drosophila read, but prior to correction, the CIGAR
            and SEQ strings should definitely match """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/drosophila_example/input_read.sam"
        genome = Fasta(f"{test_dir}/input_files/drosophila_example/chr3R.fa")
        sjDict = set()

        with open(sam, 'r') as f:
            for sam_line in f:
                if sam_line.startswith("@"):
                    continue
                else:
                    sam_line = sam_line.strip().split('\t')
                    transcript = t2.Transcript(sam_line, genome, sjDict)

        assert t2.check_seq_and_cigar_length(transcript.SEQ, transcript.CIGAR) == True

    def test_not_matching(self):
        """ This is the sam read as above but after correction with TCv2.0.1.
            The read got flagged by samtools after correction as having
            inconsistent CIGAR and SEQ fields """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/drosophila_example/bad_correction.sam"
        genome = Fasta(f"{test_dir}/input_files/drosophila_example/chr3R.fa")
        sjDict = set()

        with open(sam, 'r') as f:
            for sam_line in f:
                if sam_line.startswith("@"):
                    continue
                else:
                    sam_line = sam_line.strip().split('\t')
                    transcript = t2.Transcript(sam_line, genome, sjDict)

        assert t2.check_seq_and_cigar_length(transcript.SEQ, transcript.CIGAR) == False

    def test_perfect_match_with_introns(self):
        """ Compare SEQ and CIGAR for a transcript that is a perfect
            reference match containing introns. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_twoIntrons.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        sjDict = set()

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, sjDict)

        assert t2.check_seq_and_cigar_length(transcript.SEQ, transcript.CIGAR) == True

    def test_mismatch(self):
        """ Compare SEQ and CIGAR for a spliced transcript that contains a
            mismatch. """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/mismatch.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        sjDict = set()

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, sjDict)

        assert t2.check_seq_and_cigar_length(transcript.SEQ, transcript.CIGAR) == True


    def test_insertion_deletion_mismatch_ncsj(self):
        """ Compare SEQ and CIGAR for a transcript that contains an
           insertion, deletion, mismatch, and noncanonical splice junction in
           it. """

        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/deletion_insertion_mismatch_nc.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        sjDict = set()

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, sjDict)

        assert t2.check_seq_and_cigar_length(transcript.SEQ, transcript.CIGAR) == True
