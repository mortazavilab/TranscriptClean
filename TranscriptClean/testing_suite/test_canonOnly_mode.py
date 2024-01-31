import pytest
from pyfaidx import Fasta
import sys
import os
import subprocess
sys.path.append("..")
import TranscriptClean.transcript as t2
import TranscriptClean.TranscriptClean as TC
from TranscriptClean.dstruct import Struct
@pytest.mark.integration

class TestCanonOnly(object):
    def test_DIM_nc_not_corrected(self):
        """ This transcript has a deletion, insertion, mismatch, and NCSJ. With
            the SJ offset value of 0, the NCSJ will not be corrected, so we
            expect that in canonOnly mode, the output SAM file should be empty.
        """
        test_dir = os.path.dirname(__file__)
        # command = ["python", "../TranscriptClean.py", "--sam",
        #            f"{test_dir}/input_files/sams/deletion_insertion_mismatch_nc.sam",
        #            "--g", f"{test_dir}/input_files/hg38_chr1.fa", "-j",
        #            f"{test_dir}/input_files/GM12878_SJs_chr1.tab", "--maxLenIndel", "5",
        #            "--maxSJOffset", "0", "--correctMismatches", "True",
        #            "--correctIndels", "True", "--correctSJs", "True", "--canonOnly",
        #            "--o", f"{test_dir}/scratch/canonOnlyMode/TC"]
        command = ["transcriptclean", "--sam",
                   f"{test_dir}/input_files/sams/deletion_insertion_mismatch_nc.sam",
                   "--g", f"{test_dir}/input_files/hg38_chr1.fa", "-j",
                   f"{test_dir}/input_files/GM12878_SJs_chr1.tab", "--maxLenIndel", "5",
                   "--maxSJOffset", "0", "--correctMismatches", "True",
                   "--correctIndels", "True", "--correctSJs", "True", "--canonOnly",
                   "--o", f"{test_dir}/scratch/canonOnlyMode/TC"]

        try:
            output = subprocess.check_output(command)

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean run crashed.")


        # Now check the results

        n_lines = 0
        with open(f"{test_dir}/scratch/canonOnlyMode/TC_clean.sam", 'r') as f:
            for line in f:
                n_lines += 1

        assert n_lines == 0

        # Read logs and make sure they are OK
        expected_log = "\t".join(["c34150/f1p1/3707", "primary",
                                   "2", "0", "0",
                                   "1", "0", "0",
                                   "2", "0",
                                   "0", "1"])

        with open(f"{test_dir}/scratch/canonOnlyMode/TC_clean.log", 'r') as f:
            header = f.readline().strip()
            log = f.readline().strip()
            assert log == expected_log

        expected_TE_log = [["TranscriptID", "Position", "ErrorType", "Size",
                            "Corrected", "ReasonNotCorrected"],
                           ["c34150/f1p1/3707", "chr1_150944919",
                                 "Mismatch", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150944920",
                                 "Mismatch", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150943033_150943034",
                                 "Insertion", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150949256_150949257",
                                 "Deletion", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150950682_150950683",
                                 "Deletion", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150943994_150944918",
                                 "NC_SJ_boundary", "1", "Uncorrected",
                                 "TooFarFromAnnotJn"]]

        # Check each line of TE log
        counter = 0
        with open(f"{test_dir}/scratch/canonOnlyMode/TC_clean.TE.log", 'r') as f:
            for line in f:
                print(line)
                assert line.strip().split('\t') == expected_TE_log[counter]
                counter += 1
