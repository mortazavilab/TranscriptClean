# This script runs a series of basic tests on TranscriptClean to make sure that it is behaving as expected.
# It does not address all possible cases, but is intended to test basic functionality.

import os

def main():

    # This sam file contains a single transcript entry. The transcript has the following errors/characteristics:
    #    - a microdeletion
    #    - a microinsertion
    #    - a deletion over 5 bp
    #    - an insertion over 5 bp
    #    - a mismatch to reference that is in the same position as a SNP
    #    - a mismatch that does not overlap a known SNP
    #    - a noncanonical splice junction within 5 bp of a known junction

    # Test 1: A transcript comprised of a single exon that perfectly matches the reference genome sequence
    print "--------------------------------------------------------------------------"
    print "Running Test 1......."
    basicSamFile = "sam_files/perfectReferenceMatch_noIntrons.sam"
    runBasic(basicSamFile, "test_out/perfectReferenceMatch_noIntrons")
    print "Test 1 successful."
 
    # Test N: A transcript comprised of three exons and two introns that perfectly matches the reference genome sequence
    print "--------------------------------------------------------------------------"
    print "Running Test 2......."
    basicSamFile = "sam_files/perfectReferenceMatch_twoIntrons.sam"
    runBasic(basicSamFile, "test_out/perfectReferenceMatch_twoIntrons")
    print "Test N successful."




def runBasic(sam, outprefix):
    # Run TranscriptClean with the most basic settings: Correct microindels but nothing else.

    command = "python ../TranscriptClean.py --sam " + sam + " --genome reference_files/chr1.fa --maxLenIndel 5 --maxSJOffset 5 --outprefix " + outprefix + " > /dev/null"
    print "Running command:"
    print command

    try:
        os.system(command)
        print "Command ran successfully."
    except:
        print "ERROR: basic run failed."
 
    return 

main()
