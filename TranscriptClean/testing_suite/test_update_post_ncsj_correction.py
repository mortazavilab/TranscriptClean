import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import transcript as t2
import TranscriptClean.spliceJunction as sj
import intronBound as ib
import TranscriptClean as TC

class TestUpdatePostNCSJCorrection(object):

    def test_update(self):
        """ Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donor, acceptor, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir, chroms)
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"]
        transcript = t2.Transcript(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        # Attempt to correct the splice donor side of the junction (left)
        transcript.SEQ, transcript.CIGAR= TC.fix_one_side_of_junction(transcript.CHROM,
                                                         transcript.POS, jnNumber,
                                                         donor, 2, genome,
                                                         transcript.SEQ,
                                                         transcript.CIGAR)

        # Now test the update function
        TC.update_post_ncsj_correction(transcript, jnNumber, genome, sjDict)

        junction = transcript.spliceJunctions[jnNumber]
        assert junction.motif_code == "21"
        assert junction.isCanonical == True
        assert transcript.MD == "MD:Z:6"
        assert transcript.isCanonical == True

    def test_no_correction(self):
        """ Make sure that the attributes stay the same if no correction
            was performed
        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donor, acceptor, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir, chroms)
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"]
        transcript = t2.Transcript(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        # Now test the update function
        TC.update_post_ncsj_correction(transcript, jnNumber, genome, sjDict)

        junction = transcript.spliceJunctions[jnNumber]
        assert junction.motif_code == "0"
        assert junction.isCanonical == False
        assert transcript.MD == "MD:Z:4"
        assert transcript.isCanonical == False
