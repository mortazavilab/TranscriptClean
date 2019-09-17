import pytest
import sys
import os
import subprocess
sys.path.append("..")
import TranscriptClean as TC
@pytest.mark.unit

class TestCreateTmpSamFile(object):
    def test_create_tmp_sam(self):
        """ Create a tmp sam file from the mock header and transcripts provided,
            and make sure the correct chromosomes were detected. Then, check
            the order of the lines in the tmp file just to be sure."""

        sam_header = ["HLine1", "HLine2"]
        sam_transcripts = [ "\t".join(["read1", "mapping", "chr1", "..."]),
                            "\t".join(["read2", "mapping", "chr2", "..."]) ]
        tmp_dir = "scratch/tmp_sam_test/"                            

        fname, chroms = TC.create_tmp_sam(sam_header, sam_transcripts, tmp_dir, 
                                          process = "test")

        assert fname == "scratch/tmp_sam_test/split_uncorr_sams/test.sam"
        assert chroms == set(["chr1", "chr2"])

                    



