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
        tmp_dir = "scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

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
        tmp_dir = "scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

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
        tmp_dir = "scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

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
        tmp_dir = "scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
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


    def test_insertion(self):
        """ Compute the correct jI/jM tag for a spliced transcript that contains an
            insertion. """

        sam = "input_files/sams/insertion.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        tmp_dir = "scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, sjDict)

        correct_jM = "jM:B:c,22,22,22,22,22,22,22,22,22,22,22"
        correct_jI = ("jI:B:i,202892660,202893238,202893426,202894183,"
                      "202894283,202894590,202894750,202895521,202895718,"
                      "202896853,202896961,202909009,202909125,202911053,"
                      "202911204,202918170,202918389,202919754,202919909,"
                      "202924967,202925208,202927088")
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM


    def test_deletion_insertion_mismatch(self):
        """ Compute the correct jI/jM tag for a spliced transcript that contains an
            insertion, deletion, and mismatch. """

        sam = "input_files/sams/deletion_insertion_mismatch.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        tmp_dir = "scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, sjDict)

        correct_jM ="jM:B:c,21,21,21,21,1,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21"
        correct_jI = ("jI:B:i,154220976,154225083,154225214,154227281,"
                      "154227360,154228614,154228726,154234590,154234715,"
                      "154249375,154249438,154251040,154251319,154251480,"
                      "154251654,154253899,154254090,154254835,154254891,"
                      "154255151,154255327,154255682,154255756,154257062,"
                      "154257259,154257345,154257435,154258976,154259031,"
                      "154259947,154260030,154260891,154261110,154261591,"
                      "154261698,154266500,154266569,154268759,154268955,"
                      "154270199")
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM


    def test_insertion_deletion_mismatch_ncsj(self):
        """ Compute the correct jI/jM tag for a transcript that contains an 
           insertion, deletion, mismatch, and noncanonical splice junction in 
           it. """
 
        sam = "input_files/sams/deletion_insertion_mismatch_nc.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        tmp_dir = "scratch/test_jIjM/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        with open(sam, 'r') as f:
            sam_line = f.readline().strip()
            transcript = t2.Transcript2(sam_line, genome, sjDict)

        correct_jM ="jM:B:c,21,21,21,0,21,21,21,21,21,21,21,21,21,21,21,21,21"
        correct_jI = ("jI:B:i,150941429,150942562,150942689,150942851,150943054,"
                      "150943919,150943994,150944918,150945109,150946885,"
                      "150947013,150949121,150949279,150949366,150949526,"
                      "150950457,150951091,150951364,150951482,150959177,"
                      "150959348,150960562,150961192,150962129,150962159,"
                      "150962586,150962720,150962973,150963140,150963529,"
                      "150963742,150963994,150964084,150964246")
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM

