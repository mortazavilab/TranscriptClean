import pytest
from pyfaidx import Fasta
import sys
import os
import subprocess
sys.path.append("..")
import transcript as t2
import TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.integration

class TestCorrectTranscripts(object):
    def test_DIM_nc(self):
        """ Correct a transcript containing a deletion, insertion, mismatch,
            and noncanonical splice junction """

        # Initialize options etc.
        sam = "input_files/sams/deletion_insertion_mismatch_nc.sam"
        genome = Fasta("input_files/hg38_chr1.fa")
        sjFile = "input_files/GM12878_SJs_chr1.tab"
        tmp_dir = "scratch/example/TC_tmp/"
        os.system("mkdir -p %s" % tmp_dir)
        chroms = set(["chr1"])
        donors, acceptors, sjAnnot = TC.processSpliceAnnotation(sjFile, tmp_dir,
                                                               chroms)

        outfiles = dstruct.Struct()
        outfiles.TElog = open(tmp_dir + "DIM_nc_clean.TE.log", 'w')
        outfiles.sam = open(tmp_dir + "DIM_nc_clean.sam", 'w')
        outfiles.fasta = open(tmp_dir + "DIM_nc_clean.fasta", 'w') 
        outfiles.log = open(tmp_dir + "DIM_nc_clean.log", 'w')

        refs = dstruct.Struct()
        refs.sjAnnot = sjAnnot
        refs.genome = genome
        refs.donors = donors
        refs.acceptors = acceptors
        refs.snps = {}
        refs.deletions = {}
        refs.insertions = {}
                
        options = dstruct.Struct()
        options.maxLenIndel = 5
        options.maxSJOffset = 5
        options.correctMismatches = "true"
        options.correctIndels = "true"
        options.correctSJs = "true"
        options.primaryOnly = True
        options.canonOnly = False

        # Correct the transcript
        with open(sam, 'r') as f:
            transcripts = [f.readline().strip()]
        TC.batch_correct(transcripts, options, refs, outfiles)
           
        # Close the output files
        for handle in outfiles.values():
            handle.close() 

        # Expected transcript attributes post-correction
        correct_CIGAR = ("12M1134N126M163N202M866N74M924N191M1777N127M2109N"
                         "157M88N159M932N633M274N117M7696N170M1215N629M938N"
                         "29M428N133M254N166M390N212M253N89M163N483M")
        correct_MD = "MD:Z:3709" 
        correct_NM = "NM:i:0"
        correct_jI = ("jI:B:i,150941429,150942562,150942689,150942851,150943054,"
                      "150943919,150943994,150944917,150945109,150946885,150947013,"
                      "150949121,150949279,150949366,150949526,150950457,150951091,"
                      "150951364,150951482,150959177,150959348,150960562,150961192,"
                      "150962129,150962159,150962586,150962720,150962973,150963140,"
                      "150963529,150963742,150963994,150964084,150964246")
        correct_jM = "jM:B:c,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21"

        # Read in transcript from outfile
        with open(tmp_dir + "DIM_nc_clean.sam", 'r') as f:
            sam_line = f.readline().strip().split('\t')
        transcript = t2.Transcript(sam_line, genome, sjAnnot)
        
        assert transcript.CIGAR == correct_CIGAR
        assert transcript.MD == correct_MD
        assert transcript.NM == correct_NM 
        assert transcript.jI == correct_jI
        assert transcript.jM == correct_jM

        # Read logs and make sure they are OK
        expected_log = "\t".join(["c34150/f1p1/3707", "primary", 
                                   "2", "0", "0", 
                                   "1", "0", "0",
                                   "2", "0", 
                                   "1", "0"])
                                  
        with open(tmp_dir + "DIM_nc_clean.log", 'r') as f:
            log = f.readline().strip()
            assert log == expected_log

    def test_DIM_nc_full(self):
        """ Correct a transcript containing a deletion, insertion, mismatch,
            and noncanonical splice junction, this time by running 
            TranscriptClean from the top """

        command = ["python", "../TranscriptClean.py", "--sam",
                   "input_files/sams/deletion_insertion_mismatch_nc.sam",
                   "--g", "input_files/hg38_chr1.fa", "-j", 
                   "input_files/GM12878_SJs_chr1.tab", "--maxLenIndel", "5",
                   "--maxSJOffset", "5", "--correctMismatches", "True",
                   "--correctIndels", "True", "--correctSJs", "True", "--primaryOnly",
                   "--canonOnly", "--o", "scratch/DIM_nc_full/TC"]
                   
        try:
            output = subprocess.check_output(command)

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean run crashed.")


        # Now check the results
        # Expected transcript attributes post-correction
        correct_CIGAR = ("12M1134N126M163N202M866N74M924N191M1777N127M2109N"
                         "157M88N159M932N633M274N117M7696N170M1215N629M938N"
                         "29M428N133M254N166M390N212M253N89M163N483M")
        correct_MD = "MD:Z:3709"
        correct_NM = "NM:i:0"
        correct_jI = ("jI:B:i,150941429,150942562,150942689,150942851,150943054,"
                      "150943919,150943994,150944917,150945109,150946885,150947013,"
                      "150949121,150949279,150949366,150949526,150950457,150951091,"
                      "150951364,150951482,150959177,150959348,150960562,150961192,"
                      "150962129,150962159,150962586,150962720,150962973,150963140,"
                      "150963529,150963742,150963994,150964084,150964246")
        correct_jM = "jM:B:c,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21"

        # Read in transcript from outfile
        with open("scratch/DIM_nc_full/TC_clean.sam", 'r') as f:
            sam_line = f.readline().strip().split('\t')
       
        CIGAR = sam_line[5]
        SEQ = sam_line[9]
        MD = sam_line[14]
        NM = sam_line[13]
        jI = sam_line[16]
        jM = sam_line[15] 

        assert CIGAR == correct_CIGAR
        assert MD == correct_MD
        assert NM == correct_NM
        assert jI == correct_jI
        assert jM == correct_jM

        # Test the sequence
        with open("answer_files/DIM_nc_full/TC_clean.sam", 'r') as f:
            sam_line = f.readline().strip().split('\t')
            expected_sequence = sam_line[9]
        assert SEQ == expected_sequence

        # Read logs and make sure they are OK
        expected_log = "\t".join(["c34150/f1p1/3707", "primary",
                                   "2", "0", "0",
                                   "1", "0", "0",
                                   "2", "0",
                                   "1", "0"])

        with open("scratch/DIM_nc_full/TC_clean.log", 'r') as f:
            header = f.readline().strip()
            log = f.readline().strip()
            assert log == expected_log

        expected_TE_log = [["TranscriptID", "Position", "ErrorType", "Size",
                            "Corrected", "ReasonNotCorrected"],
                           ["c34150/f1p1/3707", "chr1_150944919",
                                 "Mismatch", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150944920",
                                 "Mismatch", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150943033_150943034",
                                 "Insertion", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150949256_150949257",
                                 "Deletion", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150950682_150950683",
                                 "Deletion", "1", "Corrected", "NA"],
                           ["c34150/f1p1/3707", "chr1_150943994_150944918",
                                 "NC_SJ_boundary", "1", "Corrected", "NA"]]

        # Check each line of TE log
        counter = 0
        with open("scratch/DIM_nc_full/TC_clean.TE.log", 'r') as f:
            for line in f:
                print(line)
                assert line.strip().split('\t') == expected_TE_log[counter]
                counter += 1

        # Make sure the report script runs without crashing
        command = ["Rscript", "../generate_report.R", 
                   "scratch/DIM_nc_full/TC"]

        try:
            output = subprocess.check_output(command)

        except Exception as e:
            print(e)
            pytest.fail("TranscriptClean report script crashed.")

        # Make sure report PDF exists and is not empty
        assert os.path.getsize("scratch/DIM_nc_full/TC_report.pdf") > 0
