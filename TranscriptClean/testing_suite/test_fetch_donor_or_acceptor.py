import pytest
from pyfaidx import Fasta
import sys
sys.path.append("..")
import spliceJunction as sj
import intronBound as ib
@pytest.mark.unit

class TestFetchDonorAcceptor(object):

    def test_fetch_donor_acceptor_plus(self):
        """ Check if SJ function correctly returns splice donor on + strand

            Toy transcript with sequence AAG|GAA, where the splice motif (GT-AG)
            is canonical.
            chr1: 23,071,357 - 23,072,126 """

        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23071360
        end = 23072123
        strand = "+"
        genome = Fasta("input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom, 
                                     start, end, strand, genome, {})

        donor = junction.get_splice_donor()
        acceptor = junction.get_splice_acceptor()

        assert donor.pos == start
        assert acceptor.pos == end  

    def test_fetch_donor_acceptor_minus(self):
        """ Check if SJ function correctly returns splice donor on - strand
        """

        transcriptID = "test_read"
        jnNumber = 0
        chrom = "chr1"
        start = 23071360
        end = 23072123
        strand = "-"
        genome = Fasta("input_files/hg38_chr1.fa")

        junction = sj.SpliceJunction(transcriptID, jnNumber, chrom,
                                     start, end, strand, genome, {})

        donor = junction.get_splice_donor()
        acceptor = junction.get_splice_acceptor()

        assert donor.pos == end
        assert acceptor.pos == start
         
