import sys
sys.path.append("..")  # noqa
import pybedtools
import subprocess
import os
import pytest
import TranscriptClean.TranscriptClean as TC
import dstruct as dstruct


@pytest.mark.integration
class TestPrepRefs(object):
    def test_genome_only(self):
        """ Make sure that the prep_refs function works under the simplest
            possible option setting: no variants or SJs provided. """
        test_dir = os.path.dirname(__file__)

        # Initialize options etc.
        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        tmp_dir = f"{test_dir}/scratch/prep_refs/genome-only/TC_tmp/"
        os.system("mkdir -p " + tmp_dir)

        options = dstruct.Struct()
        options.refGenome = f"{test_dir}/input_files/hg38_chr1.fa"
        options.tmp_dir = tmp_dir
        options.maxLenIndel = options.maxSJOffset = 5
        options.correctSJs = "false"
        options.variantFile = None
        options.sjAnnotFile = None

        header, chroms, sam_chunks = TC.split_SAM(sam, 1)
        refs = TC.prep_refs(options, sam_chunks[0], header)

        # Check that variant dicts are empty
        assert refs.snps == refs.insertions == refs.deletions == {}

        # Check that SJ bedtools and annot lookup are empty
        assert refs.donors == refs.acceptors == None
        assert refs.sjAnnot == set()

    def test_sjs(self):
        """ Genome and splice junction reference provided. Variant structs
            should still be empty. """
        test_dir = os.path.dirname(__file__)

        # Initialize options etc.
        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        tmp_dir = f"{test_dir}/scratch/prep_refs/sjs/TC_tmp/"
        os.system("mkdir -p " + tmp_dir)

        options = dstruct.Struct()
        options.refGenome = f"{test_dir}/input_files/hg38_chr1.fa"
        options.tmp_dir = tmp_dir
        options.maxLenIndel = options.maxSJOffset = 5
        options.correctSJs = "true"
        options.variantFile = None
        options.sjAnnotFile = f"{test_dir}/input_files/test_junctions.txt"

        header, chroms, sam_chunks = TC.split_SAM(sam, 1)
        refs = TC.prep_refs(options, sam_chunks[0], header)

        # Check that variant dicts are empty
        assert refs.snps == refs.insertions == refs.deletions == {}

        # Check SJ bedtools and annot lookup
        assert refs.donors.length == 3
        assert refs.acceptors.length == 2  # Same acceptor appears in 2 jns
        assert len(refs.sjAnnot) == 3

    def test_sj_corr_off(self):
        """ Splice reference provided, but correction set to off. Expected
            behavior is to skip SJ ref initialization because it would be a
            waste of time """
        test_dir = os.path.dirname(__file__)

        # Initialize options etc.
        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        tmp_dir = f"{test_dir}/scratch/prep_refs/sj_off/TC_tmp/"
        os.system("mkdir -p " + tmp_dir)

        options = dstruct.Struct()
        options.refGenome = f"{test_dir}/input_files/hg38_chr1.fa"
        options.tmp_dir = tmp_dir
        options.maxLenIndel = options.maxSJOffset = 5
        options.correctSJs = "false"
        options.variantFile = None
        options.sjAnnotFile = f"{test_dir}/input_files/test_junctions.txt"

        header, chroms, sam_chunks = TC.split_SAM(sam, 1)
        refs = TC.prep_refs(options, sam_chunks[0], header)

        # Check that variant dicts are empty
        assert refs.snps == refs.insertions == refs.deletions == {}

        # Check that SJ bedtools and annot lookup are empty
        assert refs.donors == refs.acceptors == None
        assert refs.sjAnnot == set()

    def test_variants(self):
        """ A variant file is provided """
        test_dir = os.path.dirname(__file__)

        # Initialize options etc.
        sam = f"{test_dir}/input_files/vcf_test/read_with_snps.sam"
        tmp_dir = f"{test_dir}/scratch/prep_refs/variant/TC_tmp/"
        os.system("mkdir -p " + tmp_dir)

        options = dstruct.Struct()
        options.refGenome = f"{test_dir}/input_files/hg38_chr11.fa"
        options.tmp_dir = tmp_dir
        options.maxLenIndel = options.maxSJOffset = 5
        options.correctSJs = "false"
        options.variantFile = f"{test_dir}/input_files/vcf_test/snps.vcf"
        options.sjAnnotFile = None

        header, chroms, sam_chunks = TC.split_SAM(sam, 1)
        refs = TC.prep_refs(options, sam_chunks[0], header)

        # Check that variant deletion and insertion dicts are empty
        assert len(refs.insertions) == 0
        assert len(refs.deletions) == 0
        assert len(refs.snps) > 0

        # Check that SJ bedtools and annot lookup are empty
        assert refs.donors == refs.acceptors == None
        assert refs.sjAnnot == set()
