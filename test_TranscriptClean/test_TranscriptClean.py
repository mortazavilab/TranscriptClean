# This script runs a series of basic tests on TranscriptClean to make sure that it is behaving as expected.
# It does not address all possible cases, but is intended to test basic functionality.

import os
import glob

def main():

    # This sam file contains a single transcript entry. The transcript has the following errors/characteristics:
    #    - a microdeletion
    #    - a microinsertion
    #    - a deletion over 5 bp
    #    - an insertion over 5 bp
    #    - a mismatch to reference that is in the same position as a SNP
    #    - a mismatch that does not overlap a known SNP
    #    - a noncanonical splice junction within 5 bp of a known junction

    genome = "/bio/dwyman/pacbio_f2016/data/STAR_hg38_ENCODE/hg38.fa" #"reference_files/chr1.fa"
    spliceJunctionFile = "reference_files/GM12878_SJs_chr1.tab"
    variantFile = "reference_files/GM12878_chr1.vcf"
    basicSamFile1 = "sam_files/perfectReferenceMatch_noIntrons.sam"
    basicSamFile2 = "sam_files/perfectReferenceMatch_twoIntrons.sam"
    basicSamFile1_noTags = "sam_files/perfectReferenceMatch_noIntrons_noExtraTags.sam"
    basicSamFile2_noTags = "sam_files/perfectReferenceMatch_twoIntrons_noExtraTags.sam"
    sam_DIM = "sam_files/deletion_insertion_mismatch.sam"
    sam_DIM_noTags = "sam_files/deletion_insertion_mismatch_noExtraTags.sam"
    sam_DIM_nc = "sam_files/deletion_insertion_mismatch_nc.sam"
    sam_DIM_nc_noTags = "sam_files/deletion_insertion_mismatch_nc_noExtraTags.sam"


    # A transcript comprised of a single exon (no introns) that perfectly matches the reference genome sequence
    # Correct action is to make no changes
    print "--------------------Part 1: Perfect Reference Match------------------------------------------------------"
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
  

def test_basic(sam, genome, answer, prefix):
    # Runs TranscriptClean in basic mode, then compares the output to the answer sam file using diff. If these two files are identical, the test is considered successful 
    runBasic(sam, genome, prefix)
    
    # Check whether the output matches the input
    command = "diff " + answer + " " + prefix + "_clean.sam > " + prefix + "_diff"
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

    print "\tERROR: Output is supposed to exactly match input, but this is not the case."
    return

def test_intermediate(sam, genome, sj, answer, prefix):
    # Runs TranscriptClean in intermediate mode, then compares the output to the answer sam file using diff. If these two files are identical, the test is considered successful
    runIntermediate(sam, genome, sj, prefix)

    # Check whether the output matches the input
    command = "diff " + sam + " " + prefix + "_clean.sam > " + prefix + "_diff"
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

    print "\tERROR: Output is supposed to exactly match input, but this is not the case."
    return

def test_variantAware(sam, genome, sj, variants, answer, prefix):
    # Runs TranscriptClean in variant-aware mode, then compares the output to the answer sam file using diff. If these two files are identical, the test is considered successful
    runVariantAware(sam, genome, sj, variants, prefix)

    # Check whether the output matches the input
    command = "diff " + sam + " " + prefix + "_clean.sam > " + prefix + "_diff"
    try:
        os.system(command)
        num_lines = sum(1 for line in open(prefix + "_diff"))
        if num_lines == 0:
            # Test was successful
            deleteTmpFiles(prefix)
            print "\t**Test successful**"
            return
    except: pass

    print "\tERROR: Output is supposed to exactly match input, but this is not the case."
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

        if samNew[14] == samCorrect[14]:
            print "\tMD tag: successful"
        else:
            print "\tMD tag: INCORRECT"

        if samNew[15] == samCorrect[15]:
            print "\tjM tag: successful"
        else:
            print "\tjM tag: INCORRECT"

        if samNew[16] == samCorrect[16]:
            print "\tjI tag: successful"
        else:
            print "\tjI tag: INCORRECT"

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
        print "\tERROR: Intermediate TranscriptClean run failed."
        print command
    return

deleteTmpFiles("test_out/")
print "\n"
main()
print "\n"
