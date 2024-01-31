import pytest
from pyfaidx import Fasta
import sys
import os
import subprocess
sys.path.append("..")
@pytest.mark.integration

class TestRmOpt(object):
    def test_rm_option_set(self):
        """ Check that when the --deleteTmp option is set, the TC_tmp dir is in
            fact removed. """
        test_dir = os.path.dirname(__file__)

        # Initialize options etc.
        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = f"{test_dir}/input_files/hg38_chr1.fa"
        os.system("mkdir -p scratch/option_deleteTmp")
        outprefix = f"{test_dir}/scratch/option_deleteTmp/TC"

        try:
            # subprocess.call(["python", "../TranscriptClean.py", "--sam", sam,
            #                  "--genome", genome, "--deleteTmp", "-o", outprefix])
            subprocess.call(["transcriptclean", "--sam", sam,
                             "--genome", genome, "--deleteTmp", "-o", outprefix])

            assert os.path.exists(f"{test_dir}/scratch/option_deleteTmp/TC_tmp") == False

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean crashed during test.")



    def test_rm_option_not_set(self):
        """ Check that when the --deleteTmp option is set, the TC_tmp dir is
            not removed. """
        test_dir = os.path.dirname(__file__)

        # Initialize options etc.
        sam = f"{test_dir}/input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = f"{test_dir}/input_files/hg38_chr1.fa"
        os.system("mkdir -p scratch/option_dont_deleteTmp")
        outprefix = f"{test_dir}/scratch/option_dont_deleteTmp/TC"

        try:
            subprocess.call(["transcriptclean", "--sam", sam,
                             "--genome", genome, "-o", outprefix])

            assert os.path.exists(f"{test_dir}/scratch/option_dont_deleteTmp/TC_tmp") == True

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean crashed during test.")
