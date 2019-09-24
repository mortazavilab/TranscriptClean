import pytest
from pyfasta import Fasta
import sys
import os
import subprocess
sys.path.append("..")
import transcript as t2
import TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.integration

class TestCorrectTranscripts2(object):
    def test_DIM_nc_full(self):
        """ Correct a transcript containing a variant deletion, mismatch,
            and no splice junctions by running TranscriptClean from the top.
            This samll example is derived from a larger transcript,
            m54284_181015_235905/15205058/29_3462_CCS, which originated in 
            GM12878 D9 (PB72) """

        command = ["python", "../TranscriptClean.py", "--sam",
                   "input_files/vcf_test/chr11_variantDel_andMismatch.sam",
                   "--g", "input_files/hg38_chr11.fa", "-j", 
                   "input_files/chr11_sjs.txt", 
                   "--variants", "input_files/vcf_test/chr11_variantDel.vcf",
                   "--maxLenIndel", "5",
                   "--maxSJOffset", "5", "--correctMismatches", "True",
                   "--correctIndels", "True", "--correctSJs", "True", "--primaryOnly",
                   "--o", "scratch/vDM/TC"]
                   
        try:
            output = subprocess.check_output(command)

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean run crashed.")

        # Now check the results
        # Expected transcript attributes post-correction
        correct_CIGAR = "2M1D5M"
        correct_MD = "MD:Z:2^G5"
        correct_NM = "NM:i:1"
        correct_jI = "jI:B:i,-1"
        correct_jM = "jM:B:c,-1"
        correct_seq = "GCGGGCC"

        # Read in transcript from outfile
        with open("scratch/vDM/TC_clean.sam", 'r') as f:
            header = f.readline()
            header = f.readline()
            header = f.readline()
            sam_line = f.readline().strip().split('\t')
            print(sam_line)      
 
        CIGAR = sam_line[5]
        SEQ = sam_line[9]
        MD = sam_line[-3]
        NM = sam_line[-4]
        jI = sam_line[-1]
        jM = sam_line[-2] 

        assert CIGAR == correct_CIGAR
        assert MD == correct_MD
        assert NM == correct_NM
        assert jI == correct_jI
        assert jM == correct_jM
        assert SEQ == correct_seq

        # Read logs and make sure they are OK
        expected_log = "\t".join(["m54284_181015_235905/15205058/29_3462_CCS-PB72-derived", 
                                   "primary",
                                   "0", "0", "1",
                                   "0", "0", "0",
                                   "1", "0",
                                   "0", "0"])

        with open("scratch/vDM/TC_clean.log", 'r') as f:
            header = f.readline().strip()
            log = f.readline().strip()
            assert log == expected_log

