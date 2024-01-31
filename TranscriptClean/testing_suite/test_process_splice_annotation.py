import sys
sys.path.append("..")  # noqa
import os
import TranscriptClean.TranscriptClean as TC
import pytest
import os
import subprocess
import pybedtools
import warnings


@pytest.mark.unit
class TestProcessSpliceAnnot(object):
    def test_tmp_files(self):
        """ Check that the expected tmp files are created."""

        test_dir = os.path.dirname(__file__)
        sj_file = f"{test_dir}/input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = f"{test_dir}/scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir,
                                                                chroms, process="test")

        # Check if paths of tmp files are correct
        assert os.path.exists(
            f"{test_dir}/scratch/sj_reading_test/splice_files/test_ref_splice_donors_tmp.bed")
        assert os.path.exists(
            f"{test_dir}/scratch/sj_reading_test/splice_files/test_ref_splice_acceptors_tmp.bed")
        # assert os.path.exists(
        #     f"{test_dir}/scratch/sj_reading_test/splice_files/test_ref_splice_donors_tmp.sorted.bed")
        # assert os.path.exists(
        #     f"{test_dir}/scratch/sj_reading_test/splice_files/test_ref_splice_acceptors_tmp.sorted.bed")

    def test_chrom_filtering(self):
        """ Check that only chr1 and chr2 junctions get saved"""
        test_dir = os.path.dirname(__file__)

        sj_file = f"{test_dir}/input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = f"{test_dir}/scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir,
                                                                chroms, process="test")
        # Check donor chroms
        # donor_chroms = set()
        # for pos in donor_bt:
        #     donor_chroms.add(pos.chrom)
        assert set(donor_bt.Chromosome) == chroms

        # Check acceptor chroms
        # acc_chroms = set()
        # for pos in accept_bt:
        #     acc_chroms.add(pos.chrom)
        assert set(accept_bt.Chromosome) == chroms

    def test_chrom_warning(self):
        """ Make sure the function prints a warning if no splice donors or
            acceptors are found on the provided chromosome. """
        test_dir = os.path.dirname(__file__)

        sj_file = f"{test_dir}/input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr18"])
        tmp_dir = f"{test_dir}/scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        assert pytest.warns(Warning, TC.processSpliceAnnotation, sj_file,
                            tmp_dir, chroms, process="test")

    def test_splice_donors(self):
        """ Make sure that the correct positions got labeled as splice donors """
        test_dir = os.path.dirname(__file__)

        sj_file = f"{test_dir}/input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = f"{test_dir}/scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir,
                                                                chroms, process="test")

        # Remember, file is 1-based but BedTool is 0-based
        expected_donors = set([99, 399])

        assert set(donor_bt.Start) == expected_donors

        # donors = set()
        # for donor in donor_bt:
        #     donors.add(donor.start)
        # assert donors == expected_donors

    def test_splice_acceptors(self):
        """ Make sure that the correct positions got labeled as splice acceptors """
        test_dir = os.path.dirname(__file__)

        sj_file = f"{test_dir}/input_files/toy_sjs_mixed_chroms.txt"
        chroms = set(["chr1", "chr2"])
        tmp_dir = "scratch/sj_reading_test/"
        os.system("mkdir -p " + tmp_dir)

        donor_bt, accept_bt, annot = TC.processSpliceAnnotation(sj_file, tmp_dir,
                                                                chroms, process="test")

        # Remember, file is 1-based but BedTool is 0-based
        expected_acc = set([199, 299])
        assert set(accept_bt.Start) == expected_acc

        # acceptors = set()
        # for acc in accept_bt:
        #     acceptors.add(acc.start)
        # assert acceptors == expected_acc
