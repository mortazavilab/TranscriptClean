# This script runs a series of basic tests on TranscriptClean to make sure that it is behaving as expected.
# It does not address all possible cases, but is intended to test basic functionality.

import os
import glob

def main():

    # Files used in the course of testing
    genome = "/bio/dwyman/pacbio_f2016/data/STAR_hg38_ENCODE/hg38.fa" #"reference_files/chr1.fa"
    spliceJunctionFile = "reference_files/GM12878_SJs_chr1.tab"
    variantFile = "reference_files/GM12878_chr1.vcf.gz"

    basicSamFile1 = "sam_files/perfectReferenceMatch_noIntrons.sam"
    basicSamFile2 = "sam_files/perfectReferenceMatch_twoIntrons.sam"
    basicSamFile1_noTags = "sam_files/perfectReferenceMatch_noIntrons_noExtraTags.sam"
    basicSamFile2_noTags = "sam_files/perfectReferenceMatch_twoIntrons_noExtraTags.sam"

    sam_D = "sam_files/deletion_input.sam"
    sam_D_answer = "sam_files/deletion_correctAnswer.sam"

    sam_I = "sam_files/insertion_input.sam"
    sam_I_answer = "sam_files/insertion_correctAnswer.sam"

    sam_M = "sam_files/mismatch_input.sam"
    sam_M_answer = "sam_files/mismatch_correctAnswer.sam"

    sam_DIM = "sam_files/deletion_insertion_mismatch.sam"
    sam_DIM_answer = "sam_files/deletion_insertion_mismatch_correctAnswer.sam"
    sam_DIM_noTags = "sam_files/deletion_insertion_mismatch_noExtraTags.sam"

    sam_DIM_nc = "sam_files/deletion_insertion_mismatch_nc.sam"
    sam_DIM_nc_answer = "sam_files/deletion_insertion_mismatch_nc_correctAnswer.sam"
    sam_DIM_nc_noTags = "sam_files/deletion_insertion_mismatch_nc_noExtraTags.sam"

    sam_nc_case1 = "sam_files/nc_case1.sam"
    sam_nc_case1_answer = "sam_files/nc_case1_correctAnswer.sam"
    sam_nc_case2 = "sam_files/nc_case2.sam"
    sam_nc_case2_answer = "sam_files/nc_case2_correctAnswer.sam"
    sam_nc_case3 = "sam_files/nc_case3.sam"
    sam_nc_case3_answer = "sam_files/nc_case3_correctAnswer.sam"
    sam_nc_case4 = "sam_files/nc_case4.sam"
    sam_nc_case4_answer = "sam_files/nc_case4_correctAnswer.sam"
    sam_nc_case5 = "sam_files/nc_case5.sam"
    sam_nc_case5_answer = "sam_files/nc_case5_correctAnswer.sam"

    sam_I_var = "sam_files/insertion_variant.sam"
    sam_D_var = "sam_files/deletion_variant.sam"
    sam_D_var_answer = "sam_files/deletion_variant_correctAnswer.sam"

    # A transcript comprised of a single exon (no introns) that perfectly matches the reference genome sequence
    # Correct action is to make no changes
    print "--------------------Part 1: Sanity check with perfect reference matches----------------------------------"
    print "Section A: Perfect reference match with no introns"
    print "Test input is a transcript comprised of a single exon that perfectly matches the reference genome sequence. Because the transcript does not contain any errors, the output should be identical to the original sam file in all three tests."
    print "---------------------------------------------------------------------------------------------------------"
    print "Test 1: TranscriptClean Basic Mode"
    test_basic(basicSamFile1, genome, basicSamFile1, "test_out/perfectMatch_basic_1.A.1")

    print "Test 2: TranscriptClean Intermediate Mode (Indel/Mismatch/SJ)"
    test_intermediate(basicSamFile1, genome, spliceJunctionFile, basicSamFile1, "test_out/perfectMatch_intermediate_1.A.2")

    print "Test 3: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ)"
    test_variantAware(basicSamFile1, genome, spliceJunctionFile, variantFile, basicSamFile1, "test_out/perfectMatch_variantAware_1.A.3")
 
    # A transcript comprised of three exons and two introns that perfectly matches the reference genome sequence
    # Correct action is to make no changes
    print "---------------------------------------------------------------------------------------------------------"
    print "Section B: Perfect reference match with introns"
    print "Test input is a transcript which contains two introns, but still perfectly matches the reference genome sequence. Because the transcript does not contain any errors, the output should be identical to the original sam file in all three tests."
    print "---------------------------------------------------------------------------------------------------------"
    print "Test 1: TranscriptClean Basic Mode"
    test_basic(basicSamFile2, genome, basicSamFile2, "test_out/perfectMatch_basic_1.B.1")

    print "Test 2: TranscriptClean Intermediate Mode (Indel/Mismatch/SJ)"
    test_intermediate(basicSamFile2, genome, spliceJunctionFile, basicSamFile2, "test_out/perfectMatch_intermediate_1.B.2")

    print "Test 3: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ)"
    test_variantAware(basicSamFile2, genome, spliceJunctionFile, variantFile, basicSamFile2, "test_out/perfectMatch_variantAware_1.B.3")
    print "---------------------------------------------------------------------------------------------------------"

    # Test that the code can correctly generate all extra tags (ie MD, NM, jI, jM) with and without introns
    print "--------------------Part 2: MD, NM, and jI/jM generation-------------------------------------------------"
    print "The purpose of these tests is to make sure that TranscriptClean can correctly generate all extra tags (ie MD, NM, jI, jM) for different types of transcripts."
    print "Test 1: Perfect reference match with no introns"
    test_generateTags(basicSamFile1_noTags, genome, spliceJunctionFile, basicSamFile1, "test_out/tags_2.1")

    print "Test 2: Perfect reference match with introns"
    test_generateTags(basicSamFile2_noTags, genome, spliceJunctionFile, basicSamFile2, "test_out/tags_2.2")

    print "Test 3: Transcript with introns that has an insertion, deletion, and mismatch"
    test_generateTags(sam_DIM_noTags, genome, spliceJunctionFile, sam_DIM, "test_out/tags_2.3")

    print "Test 4: Transcript with introns that has an insertion, deletion, mismatch, and a noncanonical splice junction"
    test_generateTags(sam_DIM_nc_noTags, genome, spliceJunctionFile, sam_DIM_nc, "test_out/tags_2.4")

    # Test that insertions, deletions, and mismatches are corrected properly on a limited example set.
    print "--------------------Part 3: Deletion, insertion, and mismatch correction--------------------------------"  
    print "Test 1: TranscriptClean on transcript with single deletion"
    test_basic(sam_D, genome, sam_D_answer, "test_out/deletion_basic_3.0.1")
    
    print "Test 2: TranscriptClean on transcript with single insertion"
    test_basic(sam_I, genome, sam_I_answer, "test_out/insertion_basic_3.0.2")

    print "Test 3: TranscriptClean on transcript with single mismatch"
    test_intermediate(sam_M, genome, spliceJunctionFile, sam_M_answer, "test_out/insertion_basic_3.0.3")

    print "Test 4: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ) on transcript with insertions, deletions, mismatches, and variants but canonical junctions only"
    test_variantAware(sam_DIM, genome, spliceJunctionFile, variantFile, sam_DIM_answer, "test_out/DIM_variantAware_3.0.4")

    print "Test 5: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ) on transcript with insertions, deletions, mismatches, and variants as well as a noncanonical splice junction"
    test_variantAware(sam_DIM_nc, genome, spliceJunctionFile, variantFile, sam_DIM_nc_answer, "test_out/DIM_nc_variantAware_3.0.5")

    print "Test 6: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ) on transcript with an insertion that exactly matches a variant. Correct action is to keep the transcript as-is."
    test_variantAware(sam_I_var, genome, spliceJunctionFile, variantFile, sam_I_var, "test_out/I_variantAware_3.0.6")

    print "Test 7: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ) on transcript with a deletion that exactly matches a variant. Correct action is to avoid correcting this deletion."
    test_variantAware(sam_D_var, genome, spliceJunctionFile, variantFile, sam_D_var_answer, "test_out/D_variantAware_3.0.7")

    print "--------------------Part 4: Noncanonical Splice Junction Correction-------------------------------------"
    print "Test 1: Distance from annotated junction is too large to correct"
    test_variantAware(sam_nc_case1, genome, spliceJunctionFile, variantFile, sam_nc_case1_answer, "test_out/nc_case1_4.0.1")

    print "Test 2: Exon starts late on one side and starts late on the other (dist_0 = negative, dist_1 = positive)"
    test_variantAware(sam_nc_case2, genome, spliceJunctionFile, variantFile, sam_nc_case2_answer, "test_out/nc_case2_4.0.2")

    print "Test 3: Junction is shifted in the negative direction (dist_0 and dist_1 are both negative)"
    test_variantAware(sam_nc_case3, genome, spliceJunctionFile, variantFile, sam_nc_case3_answer, "test_out/nc_case3_4.0.3")

    print "Test 4: Junction is shifted in the positive direction (dist_0 and dist_1 are both positive)"
    test_variantAware(sam_nc_case4, genome, spliceJunctionFile, variantFile, sam_nc_case4_answer, "test_out/nc_case4_4.0.4")

    print "Test 5: One side of the junction matches the annotation and the other side does not "
    test_variantAware(sam_nc_case5, genome, spliceJunctionFile, variantFile, sam_nc_case5_answer, "test_out/nc_case3_4.0.5")

def test_basic(sam, genome, answer, prefix):
    # Runs TranscriptClean in basic mode, then compares the output to the answer sam file using diff. If these two files are identical, the test is considered successful 
    runBasic(sam, genome, prefix)
    
    # Check whether the output matches the input
    command = "diff -b " + answer + " " + prefix + "_clean.sam > " + prefix + "_diff"
    try:
        os.system(command)
        num_lines = sum(1 for line in open(prefix + "_diff"))
        if num_lines == 0:
            # Test was successful
            deleteTmpFiles(prefix)
            print "\t**Test successful**"
            return 
    except:
        pass

    print "\tERROR: Output is supposed to exactly match provided file, but this is not the case."
    print command
    return

def test_intermediate(sam, genome, sj, answer, prefix):
    # Runs TranscriptClean in intermediate mode, then compares the output to the answer sam file using diff. If these two files are identical, the test is considered successful
    runIntermediate(sam, genome, sj, prefix)

    # Check whether the output matches the input
    command = "diff -b " + answer + " " + prefix + "_clean.sam > " + prefix + "_diff"
    try:
        os.system(command)
        num_lines = sum(1 for line in open(prefix + "_diff"))
        if num_lines == 0:
            # Test was successful
            deleteTmpFiles(prefix)
            print "\t**Test successful**"
            return
    except:
        pass

    print "\tERROR: Output is supposed to exactly match provided file, but this is not the case."
    print command
    return

def test_variantAware(sam, genome, sj, variants, answer, prefix):
    # Runs TranscriptClean in variant-aware mode, then compares the output to the answer sam file using diff. If these two files are identical, the test is considered successful
    runVariantAware(sam, genome, sj, variants, prefix)

    # Check whether the output matches the input
    command = "diff -b " + answer + " " + prefix + "_clean.sam > " + prefix + "_diff"
    try:
        os.system(command)
        num_lines = sum(1 for line in open(prefix + "_diff"))
        if num_lines == 0:
            # Test was successful
            deleteTmpFiles(prefix)
            print "\t**Test successful**"
            return
    except: pass

    print "\tERROR: Output is supposed to exactly match provided file, but this is not the case."
    print command
    return

def test_generateTags(sam, genome, sj, completeSam, prefix):
    # Test whether TranscriptClean generates the correct MD, NM, jI, and jM tags when these are missing from the data.
    # Correctness is checked by comparing to completeSam

    dryRun(sam, genome, sj, prefix)

    # Check that each tag matches the corresponding value in completeSam
    try:
        with open(prefix + "_clean.sam", 'r') as f:
            samNew = f.readline().strip().split("\t")
        with open(completeSam, 'r') as f:
            samCorrect = f.readline().strip().split("\t")

        if samNew[13] == samCorrect[13]:
            print "\tNM tag: successful"
        else:
            print "\tNM tag: INCORRECT"
            print "Should be: " + samCorrect[13] + ", is: " + samNew[13]

        if samNew[14] == samCorrect[14]:
            print "\tMD tag: successful"
        else:
            print "\tMD tag: INCORRECT"
            print "Should be: " + samCorrect[14] + ", is: " + samNew[14]

        if samNew[15] == samCorrect[15]:
            print "\tjM tag: successful"
        else:
            print "\tjM tag: INCORRECT"
            print "Should be: " + samCorrect[15] + ", is: " + samNew[15]

        if samNew[16] == samCorrect[16]:
            print "\tjI tag: successful"
        else:
            print "\tjI tag: INCORRECT"
            print "Should be: " + samCorrect[16] + ", is: " + samNew[16]

    except: pass   
    return

def deleteTmpFiles(prefix):
    for filename in glob.glob(prefix + "*"):
        os.remove(filename) 
    return

############################################
# Different TranscriptClean modes          #
############################################

def dryRun(sam, genome, sj, outprefix):
    # Run TranscriptClean without making corrections.
    
    command = "python ../TranscriptClean.py --sam " + sam + " --genome " + genome + " --spliceJns " + sj + "  --correctMismatches False --correctIndels False --correctSJs False --outprefix " + outprefix + " > /dev/null"
    try:
        os.system(command)
    except:
        print "\tERROR: TranscriptClean dry run failed."
        print command
    return 

def runBasic(sam, genome, outprefix):
    # Run TranscriptClean with the most basic settings: Correct microindels but nothing else.

    command = "python ../TranscriptClean.py --sam " + sam + " --genome " + genome + "  --correctMismatches False --maxLenIndel 5 --maxSJOffset 5 --outprefix " + outprefix + " > /dev/null"
    try:
        os.system(command)
    except:
        print "\tERROR: Basic TranscriptClean run failed."
        print command
    return 

def runIntermediate(sam, genome, sj, outprefix):
    # Run TranscriptClean with the intermediate settings: Correct microindels, mismatches, and noncanonical splice junctions. Correction not variant aware.

    command = "python ../TranscriptClean.py --sam " + sam + " --genome " + genome + " --spliceJns " + sj + " --correctMismatches True --maxLenIndel 5 --maxSJOffset 5 --outprefix " + outprefix + " > /dev/null"
    try:
        os.system(command)
    except:
        print "\tERROR: Intermediate TranscriptClean run failed."
        print command 
    return

def runVariantAware(sam, genome, sj, variants, outprefix):
    # Run TranscriptClean with the settings: Correct microindels, mismatches, and noncanonical splice junctions. Correction is variant aware.
    command = "python ../TranscriptClean.py --sam " + sam + " --genome " + genome + " --spliceJns " + sj + " --variants " + variants + " --correctMismatches True --maxLenIndel 5 --maxSJOffset 5 --outprefix " + outprefix + " > /dev/null"
    try:
        os.system(command)
    except:
        print "\tERROR: VariantAware TranscriptClean run failed."
        print command
    return

deleteTmpFiles("test_out/")
print "\n"
main()
print "\n"
