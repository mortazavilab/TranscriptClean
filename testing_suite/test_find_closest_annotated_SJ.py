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

    def test_find_closest_sj_plus(self):

        # Process reference junctions
        sjFile = "input_files/test_junctions.txt"
        outprefix = "scratch/test"
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, outprefix)

        # Intron bound info
        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23071350
        end = 23072124
        strand = "+"
        jnStr = "21"
        genome = Fasta("input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, jnStr, genome)

        closest_donor, closest_acceptor = TC.find_closest_ref_junction(junction, donors, acceptors)
        assert closest_donor.end == 23071360
        assert closest_acceptor.end == 23072123
        
         
