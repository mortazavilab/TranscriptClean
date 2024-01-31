import sys
import os
sys.path.append("..")  # noqa

from TranscriptClean.dstruct import Struct
import TranscriptClean.TranscriptClean as TC
import intronBound as ib
import TranscriptClean.spliceJunction as sj
import pytest
from pyfaidx import Fasta


@pytest.mark.unit
class TestFindClosestBound(object):

    def test_find_closest_splice_donor_plus(self):
        """ For a toy case with multiple donors and acceptors in close
            proximity, test whether TC can find the closest reference donor
            to the supplied intron bound.

            In this case, there is an exact match for the donor, located
            at 23071360 in 1-based coordinates and 23071359 in 0-based. """
        test_dir = os.path.dirname(__file__)

        # Process reference junctions
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23071360
        end = 23072123
        strand = "+"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, genome, sjDict)

        donor = junction.get_splice_donor()
        closest_donor = TC.find_closest_bound(donor, donors)
        assert closest_donor.start == 23071359
        assert closest_donor.end == 23071360
        assert closest_donor.dist == 0

    def test_find_closest_splice_donor_minus(self):
        """ For a toy case with multiple donors and acceptors in close
            proximity, test whether TC can find the closest reference donor
            to the supplied intron bound.

            Similar to before, there is an exact match for the donor, located
            at 23071360 in 1-based coordinates and 23071359 in 0-based."""
        test_dir = os.path.dirname(__file__)

        # Process reference junctions
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23070360
        end = 23071360
        strand = "-"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, genome, sjDict)

        donor = junction.get_splice_donor()
        closest_donor = TC.find_closest_bound(donor, donors)
        assert closest_donor.start == 23071359
        assert closest_donor.end == 23071360
        assert closest_donor.dist == 0

    def test_find_closest_splice_acceptor_plus(self):
        """ Find the closest splice acceptor, which is 17 bp upstream.
            Plus strand."""

        # Process reference junctions
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23071360
        end = 23072140
        strand = "+"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, genome, sjDict)

        acceptor = junction.get_splice_acceptor()
        closest_acceptor = TC.find_closest_bound(acceptor, acceptors)
        assert closest_acceptor.start == 23072122
        assert closest_acceptor.end == 23072123
        assert closest_acceptor.dist == -17

    def test_find_closest_splice_acceptor_minus(self):
        """ Find the closest splice acceptor, which is 1 bp downstream.
            Minus strand. Note that dist is relative to the genome, not to
            the direction of the transcript."""
        test_dir = os.path.dirname(__file__)

        # Process reference junctions
        sjFile = f"{test_dir}/input_files/test_junctions.txt"
        tmp_dir = f"{test_dir}/scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 22071331
        end = 22073331
        strand = "-"
        genome = Fasta(f"{test_dir}/input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, genome, sjDict)

        acceptor = junction.get_splice_acceptor()
        closest_acceptor = TC.find_closest_bound(acceptor, acceptors)
        assert closest_acceptor.start == 22071329
        assert closest_acceptor.end == 22071330
        assert closest_acceptor.dist == -1
