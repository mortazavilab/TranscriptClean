import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import spliceJunction as sj
import intronBound as ib
import TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.unit

class TestFindClosestSJ(object):

    def test_find_closest_splice_donor_plus(self):
        """ For a toy case with multiple donors and acceptors in close 
            proximity, test whether TC can find the closest reference donor
            to the supplied intron bound.
            
            In this case, there is an exact match for the donor, located 
            at 23071360 in 1-based coordinates and 23071359 in 0-based. """

        # Process reference junctions
        sjFile = "input_files/test_junctions.txt"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)

        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23071360
        end = 23072123
        strand = "+"
        jnStr = "21"
        genome = Fasta("input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, jnStr, genome)

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

        # Process reference junctions
        sjFile = "input_files/test_junctions.txt"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)

        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23070360
        end = 23071360
        strand = "-"
        jnStr = "25"
        genome = Fasta("input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, jnStr, genome)

        donor = junction.get_splice_donor()
        closest_donor = TC.find_closest_bound(donor, donors)
        assert closest_donor.start == 23071359
        assert closest_donor.end == 23071360
        assert closest_donor.dist == 0

    def test_find_closest_splice_acceptor_plus(self):
        """ """
        
        # Process reference junctions
        sjFile = "input_files/test_junctions.txt"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)
        
        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23071360
        end = 23072140
        strand = "+"
        jnStr = "25"
        genome = Fasta("input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, jnStr, genome)

        acceptor = junction.get_splice_acceptor()
        closest_acceptor = TC.find_closest_bound(acceptor, acceptors)
        assert closest_acceptor.start == 23072122
        assert closest_acceptor.end == 23072123
        assert closest_acceptor.dist == -17


         
        
         
