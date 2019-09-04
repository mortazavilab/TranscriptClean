import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import TranscriptClean as TC
@pytest.mark.unit

class TestComputejIjM(object):
    def test_perfect_match_no_introns(self):
        """ Compute the correct jI/jM tag for a transcript that is a perfect 
            reference match with no introns. """
        
        sam = "input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, sjDict)

        correct_jM = "jM:B:c,-1"
        correct_jI = "jI:B:i,-1"
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM

    def test_perfect_match_with_introns(self):
        """ Compute the correct jI/jM tag for a transcript that is a perfect
            reference match containing introns. """

        sam = "input_files/sams/perfectReferenceMatch_twoIntrons.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)


        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, sjDict)
        
        correct_jM = "jM:B:c,21,21"
        correct_jI = "jI:B:i,192575930,192576284,192576366,192576773"
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM

    def test_mismatch(self):
        """ Compute the correct jI/jM tag for a spliced transcript that contains a
            mismatch. """
        
        sam = "input_files/sams/mismatch.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, sjDict)

        correct_jM = "jM:B:c,21,21,21,21,21,21,21,21,21,21,21,21,21,21,23,21,21,21,21"
        correct_jI = ("jI:B:i,51730125,51745535,51745632,51746709,51746758,"
                      "51748368,51748399,51750144,51750196,51760689,51760781,"
                      "51761866,51761972,51765821,51765982,51772069,51772183,"
                      "51772604,51772724,51776832,51776919,51781163,51781336,"
                      "51782558,51782644,51783914,51784026,51784263,51784328,"
                      "51784441,51784583,51785807,51785887,51786525,51786618,"
                      "51787352,51787489,51787714")
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM


    def test_deletion(self):
        """ Compute the correct jI/jM tag for a spliced transcript that contains a
            deletion. """

        sam = "input_files/sams/deletion.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, sjDict)
    
        correct_jM = "jM:B:c,22,22,2,2,2,2,2,2,2,2"
        correct_jI = ("jI:B:i,155247956,155248049,155248304,155248404,"
                      "155248474,155250275,155250822,155251094,155251165,"
                      "155251501,155251602,155251732,155251906,155253624,"
                      "155253733,155253858,155253979,155254399,155254457,155254652")
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM

#
#    def test_insertion(self):
#        """ Compute the correct jI/jM tag for a spliced transcript that contains an
#            insertion. """
#
#        sam = "input_files/sams/insertion.sam"
#        genome = Fasta("input_files/hg38_chr1.fa")
#
#        with open(sam, 'r') as f:
#            sam_line = f.readline().strip()
#            transcript = t2.Transcript2(sam_line, genome, {})
#
#
#    def test_deletion_insertion_mismatch(self):
#        """ Compute the correct jI/jM tag for a spliced transcript that contains an
#            insertion, deletion, and mismatch. """
#
#        sam = "input_files/sams/deletion_insertion_mismatch.sam"
#        genome = Fasta("input_files/hg38_chr1.fa")
#
#        with open(sam, 'r') as f:
#            sam_line = f.readline().strip()
#            transcript = t2.Transcript2(sam_line, genome, {})
#
#
#    def test_insertion_deletion_mismatch_ncsj(self):
#        """ Compute the correct jI/jM tag for a transcript that contains an 
#           insertion, deletion, mismatch, and noncanonical splice junction in 
#           it. """
# 
#        sam = "input_files/sams/deletion_insertion_mismatch_nc.sam"
#        genome = Fasta("input_files/hg38_chr1.fa")
#
#        with open(sam, 'r') as f:
#            sam_line = f.readline().strip()
#            transcript = t2.Transcript2(sam_line, genome, {})
#
#        
