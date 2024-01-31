import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import transcriptclean.transcript as t2
import transcriptclean.TranscriptClean as TC
from transcriptclean.dstruct import Struct
@pytest.mark.unit

class TestAllSJsAnnot(object):
    def test_no_jns(self):
        """ Return transcript.allJnsAnnotated = True for a transcript without
            junctions"""
        test_dir = os.path.dirname(__file__)
        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        assert transcript.allJnsAnnotated == True
        assert transcript.isCanonical == True

    def test_two_annotated_SJs(self):
        """ Transcript with 2 junctions and each match the provided reference
        """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_twoIntrons.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        sjFile = f"{test_dir}/input_files/GM12878_SJs_chr1.tab"
        outprefix = f"{test_dir}/scratch/test"
        tmp_dir = f"{test_dir}/scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, sjDict)

        assert transcript.allJnsAnnotated == True
        assert transcript.isCanonical == True

    def test_two_annotated_SJs_without_ref(self):
        """ Same example, but no splice annot provided, so no junctions can
            show up as annotated
        """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_twoIntrons.sam"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')
            transcript = t2.Transcript(sam_line, genome, {})

        assert transcript.allJnsAnnotated == False
        assert transcript.isCanonical == True

    def test_noncanonical(self):
        """ Transcript should be noncanonical and un-annotated prior to
            correction, but be canonical and annotated afterwards """
        test_dir = os.path.dirname(__file__)

        sam = f"{test_dir}/input_files/sams/deletion_insertion_mismatch_nc.sam"
        sjFile = f"{test_dir}/input_files/GM12878_SJs_chr1.tab"
        tmp_dir = f"{test_dir}/scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        refs = Struct()
        refs.genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")
        refs.donors, refs.acceptors, refs.sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                                               chroms)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript, logInfo = TC.transcript_init(sam_line, refs.genome,
                                                     refs.sjAnnot)

        assert transcript.allJnsAnnotated == False
        assert transcript.isCanonical == False

        # Now correct the junction and retest
        upd_transcript, TE = TC.cleanNoncanonical(transcript, refs, 5, logInfo)

        assert upd_transcript.allJnsAnnotated == True
        assert upd_transcript.isCanonical == True
