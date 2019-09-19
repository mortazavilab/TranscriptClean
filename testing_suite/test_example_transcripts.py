import pytest
from pyfasta import Fasta
import sys
import os
sys.path.append("..")
import transcript as t2
import TranscriptClean as TC
import dstruct as dstruct
@pytest.mark.integration

class TestCorrectTranscripts(object):
    def test_DIM_nc(self):
        """ Correct a transcript containing a deletion, insertion, mismatch,
            and noncanonical splcie junction """

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
        outfiles.TElog = open(tmp_dir + "DIM_nc.TE.log", 'w')
        outfiles.sam = open(tmp_dir + "DIM_nc_clean.sam", 'w')
        outfiles.fasta = open(tmp_dir + "DIM_nc_clean.fasta", 'w') 
        outfiles.log = open(tmp_dir + "DIM_nc.log", 'w')

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
        options.mismatchCorrection = "true"
        options.indelCorrection = "true"
        options.sjCorrection = "true"
        options.primaryOnly = "true"

        # Correct the transcript
        with open(sam, 'r') as f:
            transcripts = [f.readline().strip()]
        TC.batch_correct(transcripts, options, refs, outfiles)
        #TC.correct_transcript(sam_line, options, refs)
           
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

