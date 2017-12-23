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

    genome = "reference_files/chr1.fa"
    spliceJunctionFile = "reference_files/GM12878_SJs_chr1.tab"
    variantFile = "reference_files/GM12878_chr1.vcf"
    basicSamFile1 = "sam_files/perfectReferenceMatch_noIntrons.sam"
    basicSamFile2 = "sam_files/perfectReferenceMatch_twoIntrons.sam"
    basicSamFile1_noTags = "sam_files/perfectReferenceMatch_noIntrons_noExtraTags.sam"
    basicSamFile2_noTags = "sam_files/perfectReferenceMatch_twoIntrons_noExtraTags.sam"

    # A transcript comprised of a single exon (no introns) that perfectly matches the reference genome sequence
    # Correct action is to make no changes
    print "--------------------Part 1: Perfect Reference Match------------------------------------------------------"
    print "Section A: Perfect reference match with no introns"
    print "---------------------------------------------------------------------------------------------------------"
    print "Test 1: TranscriptClean Basic Mode"
    test_perfectMatch_basic(basicSamFile1, genome)

    print "Test 2: TranscriptClean Intermediate Mode (Indel/Mismatch/SJ)"
    test_perfectMatch_intermediate(basicSamFile1, genome, spliceJunctionFile)

    print "Test 3: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ)"
    test_perfectMatch_variantAware(basicSamFile1, genome, spliceJunctionFile, variantFile)
 
    # A transcript comprised of three exons and two introns that perfectly matches the reference genome sequence
    # Correct action is to make no changes
    print "---------------------------------------------------------------------------------------------------------"
    print "Section B: Perfect reference match with introns"
    print "---------------------------------------------------------------------------------------------------------"
    print "Test 1: TranscriptClean Basic Mode"
    test_perfectMatch_basic(basicSamFile2, genome)

    print "Test 2: TranscriptClean Intermediate Mode (Indel/Mismatch/SJ)"
    test_perfectMatch_intermediate(basicSamFile2, genome, spliceJunctionFile)

    print "Test 3: TranscriptClean Variant-Aware Mode (Indel/Mismatch/SJ)"
    test_perfectMatch_variantAware(basicSamFile2, genome, spliceJunctionFile, variantFile)
    print "---------------------------------------------------------------------------------------------------------"

    # Test that the code can correctly generate all extra tags (ie MD, NM, jI, jM) with and without introns
    print "--------------------Part 2: MD, NM, and jI/jM generation-------------------------------------------------"
    print "Test 1: Perfect reference match with no introns"
    test_perfectMatch_generateTags(basicSamFile1_noTags, genome, basicSamFile1)

    print "Test 2: Perfect reference match with introns"
    test_perfectMatch_generateTags(basicSamFile2_noTags, genome, basicSamFile2)

def test_perfectMatch_basic(sam, genome):
    # Test input: A transcript comprised of a single exon that perfectly matches the reference genome sequence
    # Correct outcome: The output should be identical to the input sam file
 

    prefix = "test_out/perfectMatch_basic"
    runBasic(sam, genome, prefix)
    
    # Check whether the output matches the input
    command = "diff " + sam + " " + prefix + "_clean.sam > " + prefix + "_diff"
    try:
        os.system(command)
        num_lines = sum(1 for line in open(prefix + "_diff"))
        if num_lines == 0:
            # Test was successful
            deleteTmpFiles('test_out/perfectMatch_basic')
            print "\t**Test successful**"
            return 
    except:
        pass

    print "\tERROR: Output is supposed to exactly match input, but this is not the case."
    return

def test_perfectMatch_intermediate(sam, genome, sj):
    # Test input: A transcript that perfectly matches the reference genome sequence
    # Correct outcome: The output should be identical to the input sam file

    prefix = "test_out/perfectMatch_intermediate"
    runIntermediate(sam, genome, sj, prefix)

    # Check whether the output matches the input
    command = "diff " + sam + " " + prefix + "_clean.sam > " + prefix + "_diff"
    try:
        os.system(command)
        num_lines = sum(1 for line in open(prefix + "_diff"))
        if num_lines == 0:
            # Test was successful
            deleteTmpFiles('test_out/perfectMatch_intermediate')
            print "\t**Test successful**"
            return
    except:
        pass

    print "\tERROR: Output is supposed to exactly match input, but this is not the case."
    return

def test_perfectMatch_variantAware(sam, genome, sj, variants):
    # Test input: A transcript that perfectly matches the reference genome sequence
    # Correct outcome: The output should be identical to the input sam file

    prefix = "test_out/perfectMatch_variantAware"
    runVariantAware(sam, genome, sj, variants, prefix)

    # Check whether the output matches the input
    command = "diff " + sam + " " + prefix + "_clean.sam > " + prefix + "_diff"
    try:
        os.system(command)
        num_lines = sum(1 for line in open(prefix + "_diff"))
        if num_lines == 0:
            # Test was successful
            deleteTmpFiles('test_out/perfectMatch_variantAware')
            print "\t**Test successful**"
            return
    except: pass

    print "\tERROR: Output is supposed to exactly match input, but this is not the case."
    return

def test_perfectMatch_generateTags(sam, genome, completeSam):
    # Test whether TranscriptClean generates the correct MD, NM, jI, and jM tags when these are missing from the data.
    # Correctness is checked by comparing to completeSam

    prefix = "test_out/perfectMatch_tags"
    runBasic(sam, genome, prefix)

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

def runBasic(sam, genome, outprefix):
    # Run TranscriptClean with the most basic settings: Correct microindels but nothing else.

    command = "python ../TranscriptClean.py --sam " + sam + " --genome " + genome + "  --correctMismatches False --maxLenIndel 5 --maxSJOffset 5 --outprefix " + outprefix + " > /dev/null"
    try:
        os.system(command)
    except:
        print "\tERROR: Basic TranscriptClean run failed."
    return 

def runIntermediate(sam, genome, sj, outprefix):
    # Run TranscriptClean with the intermediate settings: Correct microindels, mismatches, and noncanonical splice junctions. Correction not variant aware.

    command = "python ../TranscriptClean.py --sam " + sam + " --genome " + genome + " --spliceJns " + sj + " --correctMismatches True --maxLenIndel 5 --maxSJOffset 5 --outprefix " + outprefix + " > /dev/null"
    try:
        os.system(command)
    except:
        print "\tERROR: Intermediate TranscriptClean run failed."
    return

def runVariantAware(sam, genome, sj, variants, outprefix):
    # Run TranscriptClean with the settings: Correct microindels, mismatches, and noncanonical splice junctions. Correction is variant aware.

    command = "python ../TranscriptClean.py --sam " + sam + " --genome " + genome + " --spliceJns " + sj + " --variants " + variants + " --correctMismatches True --maxLenIndel 5 --maxSJOffset 5 --outprefix " + outprefix + " > /dev/null"
    try:
        os.system(command)
    except:
        print "\tERROR: Intermediate TranscriptClean run failed."
    return

deleteTmpFiles("test_out/")
print "\n"
main()
print "\n"
