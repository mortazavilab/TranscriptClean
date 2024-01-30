import pytest
from pyfaidx import Fasta
import os
import sys
sys.path.append("..")
import transcript as t2
import spliceJunction as sj
import intronBound as ib
import TranscriptClean as TC
import dstruct as dstruct

class TestNCSJCorrection(object):

    def test_correct_ncsj(self):
        """ Toy transcript with sequence A|GAA, where the splice motif
            is noncanonical but located 2 bp from a canonical splice donor.
            chr1: 23,071,357 - 23,072,126

        """

        # Process references
        sjFile = "input_files/test_junctions.txt"
        tmp_dir = "scratch/test_ncsj/TC_tmp/"
        os.system("mkdir -p %s" % tmp_dir)
        refs = dstruct.Struct()
        chroms = set(["chr1"])
        refs.donors, refs.acceptors, refs.sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir, chroms)
        refs.genome = Fasta("input_files/hg38_chr1.fa")

        # Init transcript object
        sam_fields = ["test_read", "0", "chr1", "23071357", "255", "1M766N3M", "*",
                                "0", "0", "AGAA", "*",  "NM:i:0", "MD:Z:4"]
        transcript = t2.Transcript(sam_fields, refs.genome, refs.sjAnnot)
        jnNumber = 0
        maxDist = 5
        logInfo = TC.init_log_info(sam_fields)

        assert transcript.isCanonical == False

        # Attempt to correct the splice junction
        transcript, TE_entries = TC.cleanNoncanonical(transcript, refs, maxDist, logInfo)
        

        assert transcript.isCanonical == True
        assert transcript.spliceJunctions[jnNumber].isCanonical == True
        assert transcript.SEQ == "AAGGAA"
        assert transcript.CIGAR == "3M764N3M"
        assert transcript.MD == "MD:Z:6"
        assert logInfo.corrected_NC_SJs == 1

    def test_crash_correction(self):
        """ This is a case that is supposed to crash the NCSJ correction process,
           resulting in no correction. This is because the mapping has
           created a 7-bp micro-exon with a canonical but likely incorrect
           junction to its left, and a non-canonical junction on its right.
           Post-correction, we end up with two introns next to each other
           with a zero-length exon, which is not valid."""

        # Process references
        sjFile = "input_files/chr11_sjs.txt"
        tmp_dir = "scratch/test/TC_tmp/"
        os.system("mkdir -p %s" % tmp_dir)
        refs = dstruct.Struct()
        chroms = set(["chr11"])
        refs.donors, refs.acceptors, refs.sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir, chroms)
        refs.genome = Fasta("input_files/hg38_chr11.fa")

        sam = "input_files/sams/microexon.sam"
        with open(sam, 'r') as f:
            sam_line = f.readline().strip().split('\t')

        # Init transcript object
        transcript = t2.Transcript(sam_line, refs.genome, refs.sjAnnot)
        maxDist = 5
        logInfo = TC.init_log_info(sam_line)

        assert transcript.isCanonical == False

        # Attempt to correct the splice junction
        transcript, TE_entries = TC.cleanNoncanonical(transcript, refs, maxDist, logInfo)

        orig_CIGAR = ("1211M5612N57M464N30M2717N120M1097N23M2632N146M1225N"
                      "140M4770N72M5051N132M1513N87M567N142M3780N100M2160N"
                      "59M864N31M9891N69M1711N7M1341N47M13S")

        assert transcript.isCanonical == False
        assert transcript.MD == "MD:Z:2473"
        assert logInfo.corrected_NC_SJs == 0
        assert logInfo.uncorrected_NC_SJs == 1
        assert transcript.CIGAR == orig_CIGAR

    def test_crash_dmel(self):
        """ This is a Drosophila junction that borders a small match preceded by
            a 7 bp deletion. It is also supposed to crash correction, but did
            not in TC v2.0.1."""

        # Process references
        sjFile = "input_files/drosophila_example/chr3R_SJs.tsv"
        tmp_dir = "scratch/dmel/TC"
        os.system("mkdir -p %s" % tmp_dir)
        refs = dstruct.Struct()
        chroms = set(["chr3R"])
        refs.donors, refs.acceptors, refs.sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir, chroms)
        refs.genome = Fasta("input_files/drosophila_example/chr3R.fa")

        sam = "input_files/drosophila_example/no_SJ_corr.sam"
        with open(sam, 'r') as f:
            for sam_line in f:
                if sam_line.startswith("@"):
                    continue                
                else:
                    sam_line = sam_line.strip().split('\t')

        # Init transcript object
        transcript = t2.Transcript(sam_line, refs.genome, refs.sjAnnot)
        maxDist = 5
        logInfo = TC.init_log_info(sam_line)
        orig_CIGAR = transcript.CIGAR
        orig_seq = transcript.SEQ
        orig_MD = transcript.MD
        expected_TE = "\t".join(["m160713_133433_42182_c101000162550000001823232709161620_s1_p0/121139/11291_13013",
                                 "chr3R_14890436_14890699", "NC_SJ_boundary", "5", "Uncorrected", "Other"]) + "\n"

        assert transcript.isCanonical == False

        # Attempt to correct the splice junction
        new_transcript, TE_entries = TC.cleanNoncanonical(transcript, refs, maxDist, logInfo)

        print(TE_entries)
        assert new_transcript.isCanonical == False
        assert TE_entries == expected_TE
        assert new_transcript.MD == orig_MD
        assert logInfo.corrected_NC_SJs == 0
        assert logInfo.uncorrected_NC_SJs == 1
        assert new_transcript.CIGAR == orig_CIGAR
        assert new_transcript.SEQ == orig_seq
