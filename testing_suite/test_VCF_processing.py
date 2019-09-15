import pytest
import sys
import os
sys.path.append("..")
import TranscriptClean as TC
@pytest.mark.integration

class TestProcessVCF(object):
    def test_process_snps(self):
        """ Make sure that the processVCF function only returns the SNPs that 
            overlap the provided SAM read, and that it does not crash."""

        maxLenIndel = 5
        tmp_dir = "scratch/vcf_test"
        os.system("mkdir -p %s" % (tmp_dir))
        variant_file = "input_files/vcf_test/snps.vcf"
        sam_file = "input_files/vcf_test/read_with_snps.sam"
        snps, insertions, deletions = TC.processVCF(variant_file,
                                                    maxLenIndel,
                                                    tmp_dir,
                                                    sam_file,
                                                    add_chr = True)

        assert snps["chr11_237087"] == ["G"]
        assert "chr11_7992932" not in snps
        assert "chr18_80259190" not in snps
        assert "chr18_80259245" not in snps

