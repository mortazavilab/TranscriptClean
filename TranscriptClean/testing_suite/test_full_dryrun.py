import pytest
from pyfaidx import Fasta
import sys
import os
import subprocess
sys.path.append("..")
import transcript as t2
import TranscriptClean.TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.integration

class TestFullDryRun(object):
    def test_full(self):
        """ Run TranscriptClean in dryRun mode on a transcript that contains
            a deletion and a mismatch.
            This small example is derived from a larger transcript,
            m54284_181015_235905/15205058/29_3462_CCS, which originated in
            GM12878 D9 (PB72) """
        test_dir = os.path.dirname(__file__)

        # command = ["python", "../TranscriptClean.py", "--sam",
        #            f"{test_dir}/input_files/vcf_test/chr11_variantDel_andMismatch.sam",
        #            "--g", f"{test_dir}/input_files/hg38_chr11.fa",
        #            "--dryRun",
        #            "--o", f"{test_dir}/scratch/dryRun/"]
        command = ["transcriptclean", "--sam",
                   f"{test_dir}/input_files/vcf_test/chr11_variantDel_andMismatch.sam",
                   "--g", f"{test_dir}/input_files/hg38_chr11.fa",
                   "--dryRun",
                   "--o", f"{test_dir}/scratch/dryRun/"]

        try:
            output = subprocess.check_output(command)

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean run crashed.")

        # Read logs and make sure they are OK
        expected_log = "\t".join(["m54284_181015_235905/15205058/29_3462_CCS-PB72-derived",
                                   "primary",
                                   "NA", "1", "NA",
                                   "NA", "0", "NA",
                                   "NA", "1",
                                   "NA", "NA"])

        with open(f"{test_dir}/scratch/dryRun/TC_clean.log", 'r') as f:
            header = f.readline().strip()
            log = f.readline().strip()
            assert log == expected_log

        # Read TE log and check entries
        header = "\t".join(["TranscriptID", "Position", "ErrorType",
                "Size", "Corrected", "ReasonNotCorrected"]) + "\n"
        deletion = "\t".join(["m54284_181015_235905/15205058/29_3462_CCS-PB72-derived",
                   "chr11_207693_207694", "Deletion", "1", "Uncorrected", "DryRun"]) + "\n"
        mismatch = "\t".join(["m54284_181015_235905/15205058/29_3462_CCS-PB72-derived",
                   "chr11_207698", "Mismatch", "1", "Uncorrected", "DryRun"]) + "\n"

        expected_TE_log = [header, deletion, mismatch]

        with open(f"{test_dir}/scratch/dryRun/TC_clean.TE.log", 'r') as f:
            TE_log = f.readlines()
            assert TE_log == expected_TE_log

        # Make sure that no SAM or fasta file was created
        assert os.path.exists(f"{test_dir}/scratch/dryRun/TC_clean.sam") == False
        assert os.path.exists(f"{test_dir}/scratch/dryRun/TC_clean.fa") == False
