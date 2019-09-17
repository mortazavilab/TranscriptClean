import pytest
import sys
import os
import subprocess
import pybedtools
sys.path.append("..")
import TranscriptClean as TC
@pytest.mark.unit

class TestProcessSpliceAnnot(object):
    def test_tmp_files(self):
        """ Check that the expected tmp files are created."""

        sj_file = "input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = "scratch/sj_reading_test/"                            
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir, 
                                                                chroms, process = "test")

        # Check if paths of tmp files are correct            
        assert os.path.exists("scratch/sj_reading_test/test_ref_splice_donors_tmp.bed")
        assert os.path.exists("scratch/sj_reading_test/test_ref_splice_acceptors_tmp.bed")
        assert os.path.exists("scratch/sj_reading_test/test_ref_splice_donors_tmp.sorted.bed")
        assert os.path.exists("scratch/sj_reading_test/test_ref_splice_acceptors_tmp.sorted.bed")

    def test_chrom_filtering(self):
        """ Check that only chr1 and ch2 junctions get saved"""

        sj_file = "input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = "scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir,
                                                                chroms, process = "test")

        # Check donor chroms
        donor_chroms = set()
        for pos in donor_bt:
            donor_chroms.add(pos.chrom)
        assert donor_chroms == chroms

        # Check acceptor chroms
        acc_chroms = set()
        for pos in accept_bt:
            acc_chroms.add(pos.chrom)
        assert acc_chroms == chroms

    def test_splice_donors(self):
        """ Make sure that the correct positions got labeled as splice donors """

        sj_file = "input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = "scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir,
                                                                chroms, process = "test")

        # Remember, file is 1-based but BedTool is 0-based
        expected_donors = set([99, 399])
        donors = set()
        for donor in donor_bt:
            donors.add(donor.start)
        assert donors == expected_donors

    def test_splice_donors(self):
        """ Make sure that the correct positions got labeled as splice acceptors """

        sj_file = "input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = "scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir,
                                                                chroms, process = "test")

        # Remember, file is 1-based but BedTool is 0-based
        expected_acc = set([199, 299])
        acceptors = set()
        for acc in accept_bt:
            acceptors.add(acc.start)
        assert acceptors == expected_acc




