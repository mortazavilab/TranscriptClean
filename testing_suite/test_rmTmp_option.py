import pytest
from pyfasta import Fasta
import sys
import os
import subprocess
sys.path.append("..")
@pytest.mark.integration

class TestCorrectTranscripts(object):
    def test_rm_option_set(self):
        """ Check that when the --deleteTmp option is set, the TC_tmp dir is in
            fact removed. """

        # Initialize options etc.
        sam = "input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = "input_files/hg38_chr1.fa"
        os.system("mkdir -p scratch/option_deleteTmp")
        outprefix = "scratch/option_deleteTmp/TC"

        try:
            subprocess.call(["python", "../TranscriptClean.py", "--sam", sam,
                             "--genome", genome, "--deleteTmp", "-o", outprefix])

            assert os.path.exists("scratch/option_deleteTmp/TC_tmp") == False

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean crashed during test.")

                    

    def test_rm_option_not_set(self):
        """ Check that when the --deleteTmp option is set, the TC_tmp dir is 
            not removed. """

        # Initialize options etc.
        sam = "input_files/sams/perfectReferenceMatch_noIntrons.sam"
        genome = "input_files/hg38_chr1.fa"
        os.system("mkdir -p scratch/option_dont_deleteTmp")
        outprefix = "scratch/option_dont_deleteTmp/TC"

        try:
            subprocess.call(["python", "../TranscriptClean.py", "--sam", sam,
                             "--genome", genome, "-o", outprefix])

            assert os.path.exists("scratch/option_dont_deleteTmp/TC_tmp") == True 

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean crashed during test.")


