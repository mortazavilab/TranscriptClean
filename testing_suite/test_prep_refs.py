import pytest
import sys
import os
import subprocess
import pybedtools
sys.path.append("..")
import TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.integration

class TestPrepRefs(object):
    def test_genome_only(self):
        """ Make sure that the prep_refs function works under the simplest
            possible option setting: no variants or SJs provided. """

        # Initialize options etc.
        sam = "input_files/sams/perfectReferenceMatch_noIntrons.sam"
        tmp_dir = "scratch/prep_refs/TC_tmp/"
        os.system("mkdir -p " + tmp_dir)

        options = dstruct.Struct()
        options.refGenome = "input_files/hg38_chr1.fa"
        options.tmp_dir = tmp_dir
        options.maxLenIndel = options.maxSJOffset = 5
        options.sjCorrection = "false"
        options.variantFile = None
        options.spliceAnnot = None

        header, sam_lines = TC.split_SAM(sam)
        refs = TC.prep_refs(options, sam_lines, header) 

        # Check that variant dicts are empty
        assert refs.snps == refs.insertions == refs.deletions == {}

        # Check that SJ bedtools and annot lookup are empty
        assert refs.donors == refs.acceptors == None
        assert refs.sjDict == {}


