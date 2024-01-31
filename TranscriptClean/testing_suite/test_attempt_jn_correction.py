import pytest
from pyfaidx import Fasta
import sys
import os
sys.path.append("..")
import TranscriptClean.transcript as t2
import spliceJunction as sj
import intronBound as ib
import TranscriptClean.TranscriptClean as TC
from TranscriptClean.dstruct import Struct
@pytest.mark.unit

class TestAttemptJnCorrection(object):

    def test_too_far_away(self):
        """ A case where the NCSJ should not be corrected because it is too far
            away from the closest annotated junction relative to the maxDist
            parameter.

         Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical.
            chr1: 23,071,357 - 23,072,126
        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/input_files/test_jns/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:6"]
        transcript = t2.Transcript(sam_fields, genome, sjAnnot)
        jnNumber = 0
        maxDist = 1

        correction_status, reason, dist = TC.attempt_jn_correction(transcript,
                                                                   jnNumber,
                                                                   genome,
                                                                   donors,
                                                                   acceptors,
                                                                   sjAnnot,
                                                                   maxDist)
        assert correction_status == False
        assert reason == "TooFarFromAnnotJn"
        assert dist == 2

    def test_correct_jn(self):
        """ Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

        """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        outprefix = f"{test_dir}/input_files/test_jns/"
        tmp_dir = f"{test_dir}/input_files/test_jns/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"]
        transcript = t2.Transcript(sam_fields, genome, sjAnnot)
        jnNumber = 0
        maxDist = 5
        #donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        # Attempt to correct the splice junction
        correction_status, reason, dist = TC.attempt_jn_correction(transcript,
                                                                   jnNumber,
                                                                   genome,
                                                                   donors,
                                                                   acceptors,
                                                                   sjAnnot,
                                                                   maxDist)

        assert correction_status == True
        assert reason == "NA"
        assert dist == 2

    def test_crash(self):
        """ This is a Drosophila junction that borders a small match preceded by
            a 7 bp deletion. It is supposed to crash correction, which will result
            in a categorization of 'Other' in the log """
        test_dir = os.path.dirname(__file__)

        # Process references
        sjFile = f"{test_dir}/input_files/drosophila_example/chr3R_SJs.tsv"
        outprefix = f"{test_dir}/input_files/dmel_crash/"
        tmp_dir = f"{test_dir}/input_files/dmel_crash/TC_tmp/"
        chroms = set(["chr3R"])
        donors, acceptors, sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta(f"{test_dir}/input_files/drosophila_example/chr3R.fa")

        # Init transcript object
        sam_fields = ["test_read", "0", "chr3R", "14890420", "255", "7M7D2M264N7M", "*",
                      "0", "0", "GATCAAACAACAAGTC", "*"]
        transcript = t2.Transcript(sam_fields, genome, sjAnnot)

        jnNumber = 0
        maxDist = 5
        # Attempt to correct the splice junction
        correction_status, reason, dist = TC.attempt_jn_correction(transcript,
                                                                   jnNumber,
                                                                   genome,
                                                                   donors,
                                                                   acceptors,
                                                                   sjAnnot,
                                                                   maxDist)
        assert correction_status == False
        assert reason == "Other"
        assert dist == 5
