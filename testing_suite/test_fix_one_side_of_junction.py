import pytest
from pyfasta import Fasta
import sys
sys.path.append("..")
import transcript2 as t2
import spliceJunction as sj
import intronBound as ib
import TranscriptClean as TC

class TestFixSideOfJunction(object):

    def test_fix_donor_case1(self):
        """ Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

            So-called case # 1
        """

        # Process references
        sjFile = "input_files/test_junctions.txt"
        tmp_dir = "scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        genome = Fasta("input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = "\t".join(["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                      "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"])
        transcript = t2.Transcript2(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM, 
                                                         transcript.POS, jnNumber, 
                                                         donor, 2, genome, 
                                                         transcript.SEQ, 
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"


    def test_fix_acceptor_case2(self):
        """ Toy transcript with sequence AAG|AA, where the splice motif
            is noncanonical but located 1 bp from a canonical splice acceptor.
            chr1: 23,071,357 - 23,072,126

            So-called case #2
        """ 
        # Process references
        sjFile = "input_files/test_junctions.txt"
        tmp_dir = "scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta("input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = "\t".join(["test_read", "0", "chr1", "23071357", "255", "3M765N2M", "*",
                      "0", "0", "AAGAA", "*",  "NM:i:0", "MD:Z:5"])
        transcript = t2.Transcript2(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        acceptor = (transcript.spliceJunctions[jnNumber]).bounds[1]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM,
                                                         transcript.POS, jnNumber,
                                                         acceptor, -1, genome,
                                                         transcript.SEQ,
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"

    def test_fix_donor_case3(self):
        """ Toy transcript with sequence AAGGT|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

            So-called case #3
        """

        # Process references
        sjFile = "input_files/test_junctions.txt"
        tmp_dir = "scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta("input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = "\t".join(["test_read", "0", "chr1", "23071357", "255", "5M762N3M", "*",
                      "0", "0", "AAGGTGAA", "*",  "NM:i:0", "MD:Z:8"])
        transcript = t2.Transcript2(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        donor = (transcript.spliceJunctions[jnNumber]).bounds[0]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM,
                                                         transcript.POS, jnNumber,
                                                         donor, -2, genome,
                                                         transcript.SEQ,
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"

    def test_fix_acceptor_case4(self):
        """ Toy transcript with sequence AAG|GGAA, where the splice motif
            is noncanonical but located 1 bp from a canonical splice acceptor.
            chr1: 23,071,357 - 23,072,126

            So-called case #4
        """
        # Process references
        sjFile = "input_files/test_junctions.txt"
        tmp_dir = "scratch/test/TC_tmp/"
        chroms = set(["chr1"])
        donors, acceptors, sjDict = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)
        genome = Fasta("input_files/hg38_chr1.fa")


        # Init transcript object
        sam_fields = "\t".join(["test_read", "0", "chr1", "23071357", "255", "3M763N4M", "*",
                      "0", "0", "AAGGGAA", "*",  "NM:i:0", "MD:Z:7"])
        transcript = t2.Transcript2(sam_fields, genome, sjDict)
        jnNumber = 0
        maxDist = 5
        acceptor = (transcript.spliceJunctions[jnNumber]).bounds[1]

        # Attempt to correct the splice donor side of the junction (left)
        new_seq, new_cigar = TC.fix_one_side_of_junction(transcript.CHROM,
                                                         transcript.POS, jnNumber,
                                                         acceptor, 1, genome,
                                                         transcript.SEQ,
                                                         transcript.CIGAR)

        assert new_seq == "AAGGAA"
        assert new_cigar == "3M764N3M"

    #def test_crash_correction(self):
    #    """ This is a case that is supposed to crash the 
    #        fix_one_side_of_junction function. This is because the mapping has 
    #        created a 7-bp micro-exon with a canonical but likely incorrect 
    #        junction to its left, and a non-canonical junction on its right.
    #        Post-correction, we end up with two introns next to each other 
    #        with a zero-length exon, which is not valid."""
        #assert 1 == 0
        #sam_line = "c42860/f2p4/2485	16	chr11	2944436	255	1211M5612N57M464N30M2717N120M1097N23M2632N146M1225N140M4770N72M5051N132M1513N87M567N142M3780N100M2160N59M864N31M9891N69M1711N7M1341N47M13S	*	0	0	ctcaATGTTTCACCACTTTAATAGAATATTCTAACTATTCTGTACAAATTGAAAACACTGTATTCAGTTACAAATGTGTTCTAAGGTTAGGCCGAGCACTCTCCACAGGCCTGGTCAGTGCGGACACGGCCATCCCCGGCTGCCGGAGAGCGCCGTCACCCACTTGAAAACCCCACCCACCAGCCGCCAAGCGGTCACACCAAACCCAGGACCTTCAGACAGGATGAATGGTGGGGCCCCGCAATGGGGCTACTGAGAAAGCAGGACTTGACGCTCATACGCTCCACTGAAACGCAGGACTTCCCAGCCCAGTCCCTCAGTGGAGAAGACTGCCGAAGCCCGGCTCCGGCAGCAGGGTGGCGCCTGCGTCATGAGGACGGGCTCGCATCTTCAGCCCTGGTGGCAGGGAGCGGCGTTTCTTCCGCACAGGCCTTGCGCCTGCTAGGAAGTGGCACATCTTCCTGCTCAGGGCACCAAGGTGGTTCAGAAACGTTAAGGACGAGCCACAGCGAAAAGCCGCAGTCCTCACAGGCAAGAAGGGATAAATAAATATGAGGTGACCCGCAGCAGCTCTCACCTGGGCTGGTGTGTCACAACCCTGACCCACCCCTAAAAAAAAAAAAAAATCAAGAAGCAACATCCTAAGGAGAACAGGGCCCTACTCTACACAGCCCTTTCTGAGATGATCGGCATACAGCAGGTGATGCAGGCTGCACACTCAGCAGATTCAGCGGCTGGAAACAGCAAGTGGGTTTCTTCGGATGAAAGGGAAGAATTCAGTCCAACTGCAGGAGGGGTGGGAGAGGTTCCAGATCCTGGGAACCACATCACCAGACCTCGGCCCTTTTTGCCAAGTGACCCCCACCCCACCCTGATGTGGTCTACAGGGCCCTCCCACAGGGAAAGGCCCAGGGAAGTCCAGAGCTACAGGCACCAAGGCTGCAGAGGGTGCTGGACGAAACCTCCTATTTCTGAAATGCATTTCAGTTGCCACTGTACAAGTTAAGCAAAATAATAAGGAAAAAGGAAAAGTGAAAGTGAAAATCATGCACTTGAAAACGAGTTAGATGGAGTAAGCTCTGTCCACGGGATTGTGCTGCGGCAAGGACCGAGGCCCCGCCCACAGGCCTGGAGTCCCGACAGCCGGTCTGCCAGGCACCCGCCTCCGCTTCCTACTGCTGCTTGCATTCCGCCGGCTGGCTGGGTTCCTTCTTGGGGTTAATTTCCGCATCATCCTCGTCTTCTCCCTCCTCGTCACCTTCTAATTCCTCCTCTTCTCCTTCTTCACCTTCTTCAAAATTGTCATCATCTTCTATGGCCTCCCCAGTGAAGTACAGCACAGCCCGCGGGACTATCCGCTCACGGAAAAAGTGTCCAATTTCAAAATCAGAGGCTAATGTGAATTCAGAATCTTCATCCAGTGATTCTCCATCCCCGGATGCTTTCAATGGATTGAAGAAGTTGAAAAAGGACTCATTGGGTACTTGTTTCGTAATTGTTCTAACAGTGCCTCGACCCTTATGCTTCTGCTTTTTCTTGATGGTTTTGACAGTAACATTCTTTCCTTTCTTCCAGTCAATAGTACACCCGTCACAGTCCACAATCTCAGGACCTTCAAAGGAAAAGGGATCAGCCTTATCTGGTTCTGATTTCATCTTGTAGGTTTTTGTCAGGACTGAGTTGGTAAAGTAGTCGTTGGGTTCAAAGTGGAACTCTAACACAAAAGACATAGGCTGTCCAGGGTCAGAAAATTTCACTTTAATATCCTGCAGGTGTTTCAAGATTGGTTCATCATATTCCTGGACTAATTCACTCAGCATGTCCACATTTCTGAAGATGGTAAACCAGAACTCTGGAATTCCTTTGGGATCTGGCTCTTCAGCCGTTGCCGCTGCTTTTTCTGTGACGACTACTTTACTTTTCATGTCTCCAGCCAATTTCTCTTCCTCTTCATTTTCACTGTGCCATTCCGATTCCGCATCTGTTGGTTCAACATCGCCGGTGATAAATTCTCTTCTCTTGTCAAAGAGAGGCTGGTATAGCGCTGCATACTTTCTTTCCAAGTCATGTACCTCTTCATAGAACTTGGCTTCTATGTGAGCACATCTCACCTGAAGTTGTTTCAATGCATTAATTCTTCTTTTTACTGCTTTAGGTAAAGTTTCGATGTAGCTGGAAGGGGTGTGAGGGACATTGTCAAGTCGCTCCTGTAAAGCTGCCAGAACTCGAGGATTCTGCATCACCTGATCTGTGAGCTTTTCTGTGTTACTTGCATTTTTAGCAGCTTCCACGGAATCTGAAGGAACCCCATCTGAAAAACTGTGATCTGCCATCTGAATGTTTTTATCCCCTCAGGAGGTCCCACAAAAATTCAGCTGAGGCCAGGTAACAGCAGTTTCCTGCAGAGATCAGCTACTTCCTTTAAGACGCCTCCTgcggcagtggcggcgaCcctagcttcgctcgctttgggctgcgccgcgtgg	*	NH:i:1	HI:i:1	NM:i:0	MD:Z:2473	jM:B:c,22,22,22,22,24,22,22,22,22,22,22,22,22,22,22,0	jI:B:i,2945647,2951258,2951316,2951779,2951810,2954526,2954647,2955743,2955767,2958398,2958545,2959769,2959910,2964679,2964752,2969802,2969935,2971447,2971535,2972101,2972244,2976023,2976124,2978283,2978343,2979206,2979238,2989128,2989198,2990908,2990916,2992256"
