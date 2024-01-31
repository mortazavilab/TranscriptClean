import pytest
import sys
import os
import subprocess
sys.path.append("..")
import transcriptclean.TranscriptClean as TC
@pytest.mark.unit

class TestCreateTmpSamFile(object):
    def test_create_tmp_sam(self):
        """ Create a tmp sam file from the mock header and transcripts provided.
            Then, check the order of the lines in the tmp file just to be sure.
        """
        test_dir = os.path.dirname(__file__)

        sam_header = ["HLine1", "HLine2"]
        sam_transcripts = [ "\t".join(["read1", "mapping", "chr1", "..."]),
                            "\t".join(["read2", "mapping", "chr2", "..."]) ]
        tmp_dir = f"{test_dir}/scratch/tmp_sam_test/"

        fname, chroms = TC.create_tmp_sam(sam_header, sam_transcripts, tmp_dir,
                                          process = "test")

        assert fname == f"{test_dir}/scratch/tmp_sam_test/split_uncorr_sams/test.sam"
        assert chroms == set(["chr1", "chr2"])

        # Now check the integrity of the output file
        line_num = 0
        with open(fname, 'r') as f:
            for line in f:
                line = line.strip()
                if line_num == 0:
                    assert line == sam_header[0]
                elif line_num == 1:
                    assert line == sam_header[1]
                elif line_num == 2:
                    assert line == sam_transcripts[0]
                elif line_num == 3:
                    assert line == sam_transcripts[1]
                else:
                    pytest.fail("Output contains more lines than expected")
                line_num += 1
