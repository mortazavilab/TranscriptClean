# TranscriptClean
# Author: Dana Wyman
# 3/1/2018
# -----------------------------------------------------------------------------
# Mismatches and microindels in long reads are corrected in a  SNP-aware 
# fashion using the reference genome and a VCF file of whitelisted variants. 
# Noncanonical splice junctions can also be corrected using a file of reference 
# splice sites.

from transcript2 import Transcript2
from spliceJunction import *
from intronBound import IntronBound
from optparse import OptionParser
import pybedtools
from pyfasta import Fasta
import os
import re

def getOptions():
    parser = OptionParser()
    
    parser.add_option("--sam", "-s", dest = "sam", 
                      help = "Input SAM file containing transcripts to correct",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--genome", "-g", dest = "refGenome", 
                      help = "Reference genome fasta file. Should be the \
                      same one used during mapping to generate the sam file.",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--spliceJns", "-j", dest = "spliceAnnot", 
                      help = "Splice junction file obtained by mapping Illumina\
                       reads to the genome using STAR. More formats may be \
                      supported in the future.",
                      metavar = "FILE", type = "string", default = None)
    parser.add_option("--variants", "-v", dest = "variantFile",
                      help = "VCF formatted file of variants to avoid \
                      correcting away in the data (optional).", 
                      metavar = "FILE", type = "string", default = None )
    parser.add_option("--maxLenIndel", dest = "maxLenIndel",
                      help = "Maximum size indel to correct (Default: 5 bp)", 
                      type = "int", default = 5 )
    parser.add_option("--maxSJOffset", dest = "maxSJOffset",
                      help = "Maximum distance from annotated splice \
                      junction to correct (Default: 5 bp)", 
                      type = "int", default = 5 )
    parser.add_option("--outprefix", "-o", dest = "outprefix",
                      help = "output file prefix. '_clean' plus a file \
                      extension will be added to the end.", 
                      metavar = "FILE", type = "string", default = "out")
    parser.add_option("--correctMismatches", "-m", dest = "correctMismatches",
                      help = "If set to false, TranscriptClean will skip \
                      mismatch correction. Default: True", 
                      type = "string", default = "true" )
    parser.add_option("--correctIndels", "-i", dest = "correctIndels",
                      help = "If set to false, TranscriptClean will skip indel \
                      correction. Default: True", type = "string", 
                      default = "true")
    parser.add_option("--correctSJs", dest = "correctSJs",
                      help = "If set to false, TranscriptClean will skip \
                      splice junction correction. Default: True, but requires \
                      splice junction annotation file to work.", \
                      type = "string", default = "true") 
    parser.add_option("--dryRun", dest ="dryRun", action='store_true', 
                      help = "If this option is set, TranscriptClean will read\
                      in the sam file and record all insertions, deletions, \
                      mismatches, and noncanonical splice junctions, but it will \
                      skip correction. This is useful for checking the distribution \
                      of transcript errors in the data before running correction.")
    parser.add_option("--primaryOnly", dest ="primaryOnly", action='store_true',
                      help = "If this option is set, TranscriptClean will only \
                      output primary mappings of transcripts (ie it will filter \
                      out unmapped and multimapped lines from the SAM input.")

    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    samFile = options.sam
    genomeFile = options.refGenome
    outprefix = options.outprefix 
    variantFile = options.variantFile
    sjFile = options.spliceAnnot
    maxLenIndel = int(options.maxLenIndel)
    maxSJOffset = int(options.maxSJOffset)
    
    indelCorrection = (options.correctIndels).lower()
    mismatchCorrection = (options.correctMismatches).lower()
    sjCorrection = (options.correctSJs).lower()
    dryRun = options.dryRun
    primaryOnly = options.primaryOnly

    # Read in the reference genome. 
    print "Reading genome .............................."
    genome = Fasta(genomeFile)

    if dryRun == True:
        # Dry run mode simply catalogues all indels in the data, then exits.
        print "Dry run mode: Cataloguing indels........."
        dryRun_recordIndels(samFile, outprefix, genome)
        return

    if sjFile != None:
        print "Processing annotated splice junctions ..."
        annotatedSpliceJns, sjDict = processSpliceAnnotation(sjFile, outprefix)
    else:
        print "No splice annotation provided. Will skip splice junction correction."
        sjDict = {}

    if variantFile != None:
        print "Processing variant file ................."
        snps, insertions, deletions = processVCF(variantFile, maxLenIndel)
    else:
        print "No variant file provided. Transcript correction will not be variant-aware."
        snps = {}
        insertions = {}
        deletions = {}

    print "Processing SAM file ........................."
    oSam = open(outprefix + "_clean.sam", 'w')
    oFa = open(outprefix + "_clean.fa", 'w')
    transcriptLog = open(outprefix + "_clean.log", 'w')
    transcriptLog.write("\t".join(["TranscriptID", "Mapping", \
                        "corrected_deletions", "uncorrected_deletions", "variant_deletions", \
                        "corrected_insertions", "uncorrected_insertions", "variant_insertions", \
                        "corrected_mismatches", "variant_mismatches", \
                        "corrected_NC_SJs", "uncorrected_NC_SJs"]) + "\n")

    transcriptErrorLog = open(options.outprefix + "_clean.TE.log", 'w')
    transcriptErrorLog.write("\t".join(["TranscriptID", "Position", "ErrorType", "Size", "Corrected", "ReasonNotCorrected"]) + "\n") 

    canTranscripts, noncanTranscripts = processSAM(samFile, genome, sjDict, oSam, oFa, transcriptLog, primaryOnly) 
    if len(canTranscripts) == 0: 
        print "Note: No canonical transcripts found."
    else:
        if dryRun == True:
            dryRun_recordIndels(canTranscripts)
        else: 
            if mismatchCorrection == "true":
                print "Correcting mismatches (canonical transcripts)............"
                correctMismatches(canTranscripts, genome, snps, transcriptErrorLog)
     
            if indelCorrection == "true":
                print "Correcting insertions (canonical transcripts)............"
                correctInsertions(canTranscripts, genome, insertions, maxLenIndel, transcriptErrorLog)
                print "Correcting deletions (canonical transcripts)............"
                correctDeletions(canTranscripts, genome, deletions, maxLenIndel, transcriptErrorLog)

    if len(noncanTranscripts) == 0:
        print "Note: No noncanonical transcripts found."
    else:
        if dryRun == True:
            dryRun_recordIndels(noncanTranscripts)
        else:
            if mismatchCorrection == "true":
                print "Correcting mismatches (noncanonical transcripts)............"
                correctMismatches(noncanTranscripts, genome, snps, transcriptErrorLog)
            if indelCorrection == "true":
                print "Correcting insertions (noncanonical transcripts)............"
                correctInsertions(noncanTranscripts, genome, insertions, maxLenIndel, transcriptErrorLog)
                print "Correcting deletions (noncanonical transcripts)............"
                correctDeletions(noncanTranscripts, genome, deletions, maxLenIndel, transcriptErrorLog)
            if sjFile != None and sjCorrection == "true":
                print "Rescuing noncanonical junctions............."
                cleanNoncanonical(noncanTranscripts, annotatedSpliceJns, genome, maxSJOffset, sjDict, outprefix, transcriptErrorLog)

    print "Writing output to sam file and fasta file.................."

    # Generate the output files
    writeTranscriptOutput(canTranscripts, sjDict, oSam, oFa, transcriptLog, genome)
    writeTranscriptOutput(noncanTranscripts, sjDict, oSam, oFa, transcriptLog, genome)

    oSam.close()
    oFa.close()
    transcriptLog.close()
    transcriptErrorLog.close()

def writeTranscriptOutput(transcripts, spliceAnnot, outSam, outFa, transcriptLog, genome):
    """Given a dict containing transcripts, write them to sam and fasta output files. Also write
       a log file that tracks every transcript. """

    for t in transcripts.keys():
        currTranscript = transcripts[t]
        outSam.write(Transcript2.printableSAM(currTranscript, genome, spliceAnnot) + "\n")
        outFa.write(Transcript2.printableFa(currTranscript) + "\n")
        transcriptLog.write("\t".join([str(i) for i in currTranscript.logInfo]) + "\n")
    return

def processSAM(sam, genome, spliceAnnot, outSam, outFa, outTLog, primaryOnly):
    """ Extracts the SAM header (for use in the output files at the end) and 
    creates a Transcript object for every transcript in the SAM file. 
    Transcripts are returned two separate dicts: one for canonical transcripts 
    and one for noncanonical. Some transcripts are unmapped or are non-primary 
    multimapping alignments. These are outputted directly and are not corrected."""

    canTranscripts = {}
    noncanTranscripts = {}

    with open(sam, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith("@"): # header line
                outSam.write(line + "\n")
                continue
            
            t = Transcript2(line, genome, spliceAnnot)

            # Unmapped/multimapper cases
            if t.mapping != 1:
                if primaryOnly == True:
                    continue
                if t.mapping == 0:
                    logInfo = [t.QNAME, "unmapped", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
                    outFa.write(Transcript2.printableFa(t) + "\n")
                elif t.mapping == 2:
                    logInfo = [t.QNAME, "non-primary", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]               
                outSam.write(Transcript2.printableSAM(t, genome, spliceAnnot) + "\n")
                outTLog.write("\t".join(logInfo) + "\n") 
                continue

            # Assign to canonical or noncanonical group based on splice jns
            if t.isCanonical == True:
                canTranscripts[t.QNAME] = t
            elif t.isCanonical == False:
                noncanTranscripts[t.QNAME] = t

    return canTranscripts, noncanTranscripts 

def processSpliceAnnotation(annotFile, outprefix):
    """ Reads in the tab-separated STAR splice junction file and creates a 
        bedtools object. Also creates a dict (annot) to allow easy lookup 
        to find out if a splice junction is annotated or not """

    bedstr = ""
    annot = {}
    fName = outprefix + "_ref_splice_jns_tmp.bed"
    fNameSorted = outprefix + "_ref_splice_jns_tmp.sorted.bed"
    o = open(fName, 'w')
    with open(annotFile, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            strand = "."
            if fields[3] == "1": strand = "+"
            if fields[3] == "2": strand = "-"

            intronMotif = int(fields[4])
            annotated = int(fields[5])
            uniqueReads = fields[6]
            multiReads = fields[7]
            maxOverhang = fields[8]

            # Process only canonical splice junctions or annotated noncanonical junctions
            if (intronMotif != 0) or (intronMotif == 0 and annotated == 1):
                # Make one bed entry for each end of the junction
                bed1 = "\t".join([chrom, str(start - 1), str(start), ".", uniqueReads, strand])
                bed2 = "\t".join([chrom, str(end - 1), str(end), ".", uniqueReads, strand])
                o.write(bed1 + "\n")
                o.write(bed2 + "\n")

                annot["_".join([chrom, str(start)])] = 1
                annot["_".join([chrom, str(end)])] = 1
    o.close()

    # Convert bed file into BedTool object
    os.system('bedtools sort -i ' + fName + ' >' + fNameSorted)
    spliceJnBedTool = pybedtools.BedTool(fNameSorted)
    os.system("rm " + fName)
    return spliceJnBedTool, annot

def processVCF(vcf, maxLen):
    """ This function reads in variants from a VCF file and stores them. SNPs 
        are stored in a dictionary by position. Indels are stored in their own 
        dictionary by start and end position."""
    SNPs = {}
    insertions = {}
    deletions = {}

    # Check if SNP file is gzipped. If it is, access contents with zcat
    if vcf.endswith(".gz"):
        tmpFile = "tmp.snps"
        os.system("zcat " + vcf + " > " +  tmpFile)
        vcf = tmpFile

    with open(vcf, 'r') as f:
        for line in f:
            line = line.strip()
    
            if line.startswith("#"): continue
            fields = line.split()
            
            chrom = fields[0]
            if not chrom.startswith("chr"): chrom = "chr" + chrom
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(",")
            refLen = len(ref)

            # It is possible for a single position to have more than one 
            # type of alternate allele (ie mismatch or indel). So it is 
            # necessary to treat each alternate allele separately.
            for allele in alt:
                altLen = len(allele)

                # SNP case
                if refLen == altLen == 1:
                    ID = chrom + "_" + str(pos)
                    if ID not in SNPs:
                        SNPs[ID] = [allele]
                    else:
                        SNPs[ID].append(allele)
                # Insertion/Deletion
                else:
                    size = abs(refLen - altLen)
                    # Only store indels of correctable size
                    if size > maxLen: continue 
                    # Positions in VCF files are one-based. 
                    # Inserton/deletion sequences always start with the 
                    # preceding normal reference base, so we need to add 
                    # one to pos in order to match with transcript positions.
                    # This principle holds for both insertions and deletions.
                    actPos = pos + 1
                    allele = allele[1:]
                    ID = "_".join([ chrom, str(actPos), str(actPos + size - 1)])

                    if refLen - altLen < 0: # Insertion
                        if ID not in insertions:
                            insertions[ID] = [allele]
                        else:
                            insertions[ID].append(allele)
                    elif refLen - altLen > 0: # Deletion
                        deletions[ID] = 1

    return SNPs, insertions, deletions

def correctInsertions(transcripts, genome, variants, maxLen, transcriptErrorLog):
    """ Corrects insertions up to size maxLen using the reference genome. 
        If a variant file was provided, correction will be SNP-aware. """

    for t in transcripts.keys():
        t = transcripts[t]
  
        # Only correct primary mapping transcripts
        if t.mapping != 1: continue

        origSeq = t.SEQ
        origCIGAR = t.CIGAR
         
        cigarOps,cigarCounts = t.splitCIGAR() 

        # Check for insertions. If none are present, we can skip this transcript
        if "I" not in origCIGAR : continue

        newCIGAR = ""
        newSeq = ""
        MVal = 0
        seqPos = 0

        # Start at position in the genome where the transcript starts.
        genomePos = t.POS 
    
        # Iterate over operations to sequence and repair insertions
        for op,ct in zip(cigarOps, cigarCounts):
 
            currPos = t.CHROM + ":" + str(genomePos) + "-" + str(genomePos + ct - 1)          
  
            if op == "M":
                 newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                 MVal += ct
                 seqPos += ct
                 genomePos += ct
             
            if op == "I":
                ID = "_".join([t.CHROM, str(genomePos), str(genomePos + ct - 1)])

                # Only insertions of a given size are corrected
                if ct <= maxLen:
                    # Check if the insertion is in the optional variant catalog.
                    if ID in variants:
                        # The insertion perfectly matches a variant position. 
                        # Leave the sequence alone if it matches an allele sequence.
                        currSeq = origSeq[seqPos:seqPos + ct]
                        if currSeq in variants[ID]:
                            Transcript2.addVariantInsertion(t)
                            errorEntry = "\t".join([t.QNAME, ID, "Insertion", str(ct), "Uncorrected", "VariantMatch"])
                            transcriptErrorLog.write(errorEntry + "\n")  
 
                            # Leave insertion in 
                            MVal, newCIGAR = endMatch(MVal, newCIGAR)
                            newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                            newCIGAR = newCIGAR + str(ct) + op
                            seqPos += ct
                            continue

                    # Correct insertion
                    errorEntry = "\t".join([t.QNAME, ID, "Insertion", str(ct), "Corrected", "NA"])
                    Transcript2.addCorrectedInsertion(t)
                    transcriptErrorLog.write(errorEntry + "\n")
                    # Subtract the inserted bases by skipping them. 
                    # GenomePos stays the same, as does MVal
                    seqPos += ct
                else: # Move on without correcting insertion because it is too big
                    errorEntry = "\t".join([t.QNAME, ID, "Insertion", str(ct), "Uncorrected", "TooLarge"])
                    Transcript2.addUncorrectedInsertion(t)
                    transcriptErrorLog.write(errorEntry + "\n")
                    MVal, newCIGAR = endMatch(MVal, newCIGAR)
                    newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                    newCIGAR = newCIGAR + str(ct) + op
                    seqPos += ct

            if op == "S":
                # End any ongoing match
                MVal, newCIGAR = endMatch(MVal, newCIGAR)
                newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                newCIGAR = newCIGAR + str(ct) + op
                seqPos += ct

            # N, H, and D operations are cases where the transcript sequence 
            # is missing bases that are in the reference genome.
            if op in ["N", "H", "D"]:
                # End any ongoing match
                MVal, newCIGAR = endMatch(MVal, newCIGAR)
                genomePos += ct
                newCIGAR = newCIGAR + str(ct) + op

        # End any ongoing match
        MVal, newCIGAR = endMatch(MVal, newCIGAR)

        # Update transcript 
        t.CIGAR = newCIGAR
        t.SEQ = newSeq
        
    return    

def correctDeletions(transcripts, genome, variants, maxLen, transcriptErrorLog):
    """ Corrects deletions up to size maxLen using the reference genome. 
        If a variant file was provided, correction will be SNP-aware."""

    for t in transcripts.keys():
        t = transcripts[t]

        origSeq = t.SEQ
        origCIGAR = t.CIGAR

        cigarOps,cigarCounts = t.splitCIGAR()

        # Check for deletions. If none are present, we can skip this transcript
        if "D" not in origCIGAR : continue

        newCIGAR = ""
        newSeq = ""
        MVal = 0
        seqPos = 0

        # Start at position in the genome where the transcript starts.
        genomePos = t.POS

        # Iterate over operations to sequence and repair mismatches and microindels
        for op,ct in zip(cigarOps, cigarCounts):
 
            currPos = t.CHROM + ":" + str(genomePos) + "-" + str(genomePos + ct - 1) 

            if op == "M":
                 newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                 MVal += ct
                 seqPos += ct
                 genomePos += ct

            if op == "D":
                ID = "_".join([t.CHROM, str(genomePos), str(genomePos + ct - 1)])
                if ct <= maxLen:

                    # Check if the deletion is in the optional variant catalog.
                    if ID in variants:
                        # The deletion perfectly matches a deletion variant. Leave the deletion in.
                        currSeq = genome.sequence({'chr': t.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                        errorEntry = "\t".join([t.QNAME, ID, "Deletion", str(ct), "Uncorrected", "VariantMatch"])
                        Transcript2.addVariantDeletion(t)
                        transcriptErrorLog.write(errorEntry + "\n")

                        MVal, newCIGAR = endMatch(MVal, newCIGAR)
                        genomePos += ct
                        newCIGAR = newCIGAR + str(ct) + op
                        continue

                    # Correct deletion if we're not in variant-aware mode
                    errorEntry = "\t".join([t.QNAME, ID, "Deletion", str(ct), "Corrected", "NA"])
                    Transcript2.addCorrectedDeletion(t)
                    transcriptErrorLog.write(errorEntry + "\n")

                    # Add the missing reference bases
                    refBases = genome.sequence({'chr': t.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                    newSeq = newSeq + refBases
                    genomePos += ct
                    MVal += ct

                # Deletion is too big to fix
                else:
                    errorEntry = "\t".join([t.QNAME, ID, "Deletion", str(ct), "Uncorrected", "TooLarge"])
                    Transcript2.addUncorrectedDeletion(t)
                    transcriptErrorLog.write(errorEntry + "\n")

                    # End any ongoing match
                    MVal, newCIGAR = endMatch(MVal, newCIGAR)
                    newCIGAR = newCIGAR + str(ct) + op
                    genomePos += ct

            # S and I operations are cases where the transcript sequence contains bases that are not in the reference genome. 
            if op in ["S", "I"]:
                # End any ongoing match
                MVal, newCIGAR = endMatch(MVal, newCIGAR)
                newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                newCIGAR = newCIGAR + str(ct) + op
                seqPos += ct

            # N and H operations are cases where the transcript sequence is missing bases that are in the reference genome. 
            if op in ["N", "H"]:
                # End any ongoing match
                MVal, newCIGAR = endMatch(MVal, newCIGAR)
                genomePos += ct
                newCIGAR = newCIGAR + str(ct) + op

        # End any ongoing match
        MVal, newCIGAR = endMatch(MVal, newCIGAR)

        t.CIGAR = newCIGAR
        t.SEQ = newSeq
        t.NM, t.MD = t.getNMandMDFlags(genome)
    return 


def correctMismatches(transcripts, genome, variants, transcriptErrorLog):
    """ This function corrects mismatches in the sequences. If a variant file was
        provided, correction will be SNP-aware."""
   
    for transcript in transcripts.keys():
        transcript = transcripts[transcript]

        origSeq = transcript.SEQ
        origCIGAR = transcript.CIGAR
        origMD = transcript.MD

        # Check for mismatches. If none are present, we can skip this transcript
        if any(i in origMD.upper() for i in 'ACTGN') == False : continue
        newCIGAR = ""
        newSeq = ""
        MVal = 0
        seqPos = 0
        genomePos = transcript.POS

        # Merge CIGAR and MD tag information so that we have the locations of all insertions, deletions, and mismatches
        mergeOperations, mergeCounts = transcript.mergeMDwithCIGAR()

        # Iterate over operations to sequence and repair mismatches and microindels
        for op,ct in zip(mergeOperations, mergeCounts):
            if op == "M":
                 newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                 MVal += ct
                 seqPos += ct
                 genomePos += ct

            # This denotes a mismatch
            if op == "X":
                # Check if the position matches a variant in the optional variant catalog. 
                ID = transcript.CHROM + "_" + str(genomePos)
                if ID in variants:
                # Mismatch position matches a known variant. If the base matches an alternative allele, do not correct it.
                    currBase = origSeq[seqPos]
                    if currBase in variants[ID]: 
                        errorEntry = "\t".join([transcript.QNAME, ID, "Mismatch", str(ct), "Uncorrected", "VariantMatch"])
                        transcriptErrorLog.write(errorEntry + "\n")
                        Transcript2.addVariantMismatch(transcript)

                        # Keep the base as-is
                        newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                        MVal += ct
                        seqPos += ct
                        genomePos += ct
                        continue
                # Otherwise, correct the mismatch to reference base
                errorEntry = "\t".join([transcript.QNAME, ID, "Mismatch", str(ct), "Corrected", "NA"])
                transcriptErrorLog.write(errorEntry + "\n")
                Transcript2.addCorrectedMismatch(transcript)

                # Change sequence base to the reference base at this position
                newSeq = newSeq + genome.sequence({'chr': transcript.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                seqPos += ct # skip the original sequence base
                genomePos += ct # advance the genome position
                MVal += ct

            if op in ["S", "I"]:
                # End any ongoing match
                MVal, newCIGAR = endMatch(MVal, newCIGAR)
                newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                newCIGAR = newCIGAR + str(ct) + op
                seqPos += ct

            if op in ["N", "H", "D"]:
                # End any ongoing match
                MVal, newCIGAR = endMatch(MVal, newCIGAR)
                genomePos += ct
                newCIGAR = newCIGAR + str(ct) + op

        # End any ongoing match
        MVal, newCIGAR = endMatch(MVal, newCIGAR)

        transcript.CIGAR = newCIGAR
        transcript.SEQ = newSeq
        transcript.NM, transcript.MD = transcript.getNMandMDFlags(genome)
    return 

                        
def endMatch(MVal, newCIGAR):
    """ Function to end an ongoing match during error correction, updating new 
        CIGAR and MD tag strings as needed. MVal keeps track of how many 
        matched bases we have at the moment so that when errors are fixed, 
        adjoining matches can be combined. When an intron or unfixable error 
        is encountered, it is necessary to end the match by outputting the 
        current MVal and then setting the variable to zero."""

    if MVal > 0:
        newCIGAR = newCIGAR + str(MVal) + "M"
        MVal = 0

    return MVal, newCIGAR

def cleanNoncanonical(transcripts, annotatedJunctions, genome, n, spliceAnnot, outprefix, transcriptErrorLog):
    """ Iterate over noncanonical transcripts. Determine whether each end is 
        within n basepairs of an annotated junction. If it is, run the rescue 
        function on it. If not, discard the transcript."""

    fName = outprefix + "_noncanonicalJnHalves_tmp.bed"
    fNameSorted = outprefix + "_sorted_noncanonicalJnHalves_tmp.bed"
    o = open(fName, 'w')
    for tID in transcripts.keys():
        t = transcripts[tID]

        bounds = Transcript2.getAllIntronBounds(t)

        for b in bounds:
            if b.isCanonical == True:
                continue

            # Get BedTool object for start of junction
            pos = IntronBound.getBED(b)
            o.write(pos + "\n")

    o.close()
    os.system('bedtools sort -i ' + fName + ' > ' + fNameSorted)
    nc = pybedtools.BedTool(fNameSorted)
    jnMatchFile = outprefix + "_jnMatches_tmp.bed"
    (nc.closest(annotatedJunctions, s=True, D="ref", t="first")).saveas(jnMatchFile)
    jnMatchFileSorted = outprefix + "_jnMatches_sorted_tmp.bed"
    os.system("sort -k 4,4 " + jnMatchFile  + " > " + jnMatchFileSorted)
    
    # Cleanup 
    os.system("rm " + fName)
    os.system("rm " + fNameSorted)

    # Iterate over noncanonical splice junction boundaries and their closest
    # canonical match. We process two consecutive lines at once because these 
    # represent each side of the jn.
    with open(jnMatchFileSorted, 'r') as f:
        for junction_half_0 in f:
            junction_half_1 = f.next()
            
            junction_half_0 = junction_half_0.strip().split("\t")
            junction_half_1 = junction_half_1.strip().split("\t")
            
            # Distances from annotated junction
            dist_0 = int(junction_half_0[-1])
            dist_1 = int(junction_half_1[-1])

            # Get information about the transcript that the junction belongs to
            transcriptID, spliceJnNum = junction_half_0[3].split("__")[:2]
            currTranscript = transcripts[transcriptID]
            currJunction = currTranscript.spliceJunctions[int(spliceJnNum)]
            ID = "_".join([currJunction.chrom, str(currJunction.start), str(currJunction.end)]) 
            
            # Only attempt to rescue junction boundaries that are within n bp of an annotated junction
            combinedDist = combinedJunctionDist(dist_0, dist_1)
            if combinedDist > n or abs(dist_0) > 2*n or abs(dist_1) > 2*n:
                errorEntry = "\t".join([currTranscript.QNAME, ID, "NC_SJ_boundary", str(combinedDist), "Uncorrected", "TooFarFromAnnotJn"])
                transcriptErrorLog.write(errorEntry + "\n")
                Transcript2.addUncorrected_NC_SJ(currTranscript)
            else:
                for jn in [junction_half_0, junction_half_1]: # Perform correction 
                    side = jn[3].split("__")[-1]
                    currIntronBound = currJunction.bounds[int(side)]
                    currDist = int(jn[-1])
                    
                    # It is possible for one side of the junction to match the annotation. If so, do not change this side
                    if currDist == 0:
                        currIntronBound.isCanonical = True
                        SpliceJunction.recheckPosition(currJunction)
                        SpliceJunction.recheckJnStr(currJunction, genome, spliceAnnot)       
                    else:
                        rescueNoncanonicalJunction(currTranscript, currJunction, currIntronBound, currDist, genome, spliceAnnot)

                errorEntry = "\t".join([currTranscript.QNAME, ID, "NC_SJ_boundary", str(combinedDist), "Corrected", "NA"])
                transcriptErrorLog.write(errorEntry + "\n")
                Transcript2.addCorrected_NC_SJ(currTranscript)

                currTranscript.NM, currTranscript.MD = currTranscript.getNMandMDFlags(genome)
    return

def combinedJunctionDist(dist_0, dist_1):
    """Computes the combined genomic distance of two splice junction ends
    from the closest annotated junctions. In essence, it is finding the 
    size indel that could have created the discrepancy between the 
    reference and transcript junctions.    

    Examples ('|' character respresents end of exon):

        Reference:     ----->|          |<-----
        Transcript:      ----->|      |<-----
            dist_0 = -2, dist_1 = +2, combined dist = 4

        Reference:     ----->|          |<-----
        Transcript:    ----->|        |<-----
            dist_0 = 0, dist_1 = +2, combined dist = 2

        Reference:     ----->|          |<-----
        Transcript:   ----->|       |<-----
            dist_0 = +1, dist_1 = +4, combined dist = 3

    """
    # If dist_0 and dist_1 have different signs, the combined distance is
    # the sum of their absolute values
    if dist_0*dist_1 <= 0:
        combined_dist = abs(dist_0) + abs(dist_1)
    else:
        combined_dist = abs(abs(dist_0) - abs(dist_1))
    return combined_dist


def rescueNoncanonicalJunction(transcript, spliceJn, intronBound, d, genome, spliceAnnot):
    """ Corrects a noncanonical splice junction, including the sequence and CIGAR string.
        Rechecks junction motif and updates the transcript. NM and MD are not updated-
        this happens in cleanNoncanonical, the function that calls this one."""
    oldCIGAR = transcript.CIGAR
    seq = transcript.SEQ

    currExonStr = ""
    currExonCIGAR = ""
    currSeqIndex = 0
    operations, counts = splitCIGAR(transcript.CIGAR)
    exonSeqs = []
    exonCIGARs = []
    intronCIGARs = []

    # First, use the old CIGAR string to segment the sequence by exon
    # Also, group the CIGAR string by exon
    for op, ct in zip(operations, counts):
        if op == "N":
            exonSeqs.append(currExonStr)
            intronCIGARs.append(ct)
            exonCIGARs.append(currExonCIGAR)
            currExonStr = ""
            currExonCIGAR = ""
        elif op in [ "M", "S", "I"]:
            currExonStr = currExonStr + str(seq[currSeqIndex:currSeqIndex+ct])
            currExonCIGAR = currExonCIGAR + str(ct) + op
            currSeqIndex += ct
        else:
            currExonCIGAR = currExonCIGAR + str(ct) + op
    exonSeqs.append(currExonStr)
    exonCIGARs.append(currExonCIGAR)

    # Now go and fix the specific splice junction piece
    targetJn = spliceJn.jnNumber

    if intronBound.bound == 0:
        targetExon = targetJn
        exon = exonSeqs[targetExon]
        if d > 0: # Need to add d bases from reference to end of exon. Case 1.
            # For CIGAR string, 
            exonEnd = intronBound.pos - 1
            seqIndex = exonEnd - transcript.POS + 1
            refAdd = genome.sequence({'chr': transcript.CHROM, 'start': exonEnd + 1, 'stop': exonEnd + d}, one_based=True)
            exonSeqs[targetExon] = exon + refAdd
            intronBound.pos += d
        if d < 0: # Need to subtract from end of exon sequence. Case 3
            exonSeqs[targetExon] = exon[0:d]
            intronBound.pos += d
        intronCIGARs[targetJn] -= d
        exonCIGARs[targetExon] = editExonCIGAR(exonCIGARs[targetExon], -1, d)
    else:
        targetExon = targetJn + 1
        exon = exonSeqs[targetExon]
        if d < 0: # Need to add d bases from reference to start of exon sequence. Case 2.
            exonStart = intronBound.pos + 1
            seqIndex = exonStart - transcript.POS + 1
            refAdd = genome.sequence({'chr': transcript.CHROM, 'start': exonStart - abs(d), 'stop': exonStart - 1}, one_based=True)
            exonSeqs[targetExon] = refAdd + exon
            intronBound.pos += d
        if d > 0: # Need to subtract from start of exon sequence. Case 4
            exonSeqs[targetExon] = exon[d:]
            intronBound.pos += d

        # Modify exon string
        exonCIGARs[targetExon] = editExonCIGAR(exonCIGARs[targetExon], 0, -d)
        intronCIGARs[targetJn] += d
    transcript.SEQ = str(''.join(exonSeqs))

    # Paste together the new CIGAR string
    newCIGAR = ""
    for i in range(0,len(intronCIGARs)):
        newCIGAR = newCIGAR + exonCIGARs[i] + str(intronCIGARs[i]) + "N"
    newCIGAR = newCIGAR + exonCIGARs[-1]
    transcript.CIGAR = newCIGAR    
    intronBound.isCanonical = True
    SpliceJunction.recheckPosition(spliceJn)
    SpliceJunction.recheckJnStr(spliceJn, genome, spliceAnnot) 

    return

 
def splitCIGAR(CIGAR):
    """ Takes CIGAR string from SAM and splits it into two lists: one with 
        capital letters (match operators), and one with the number of bases """

    matchTypes = re.sub('[0-9]', " ", CIGAR).split()
    matchCounts = re.sub('[A-Z]', " ", CIGAR).split()
    matchCounts = [int(i) for i in matchCounts]

    return matchTypes, matchCounts

def editExonCIGAR(exon, side, nBases):
    """ Given an exon CIGAR string, this function adds or subtracts bases from
        the start (side = 0) or end (side = -1)"""

    operations, counts = splitCIGAR(exon)
    if side == -1:
        operations = operations[::-1]
        counts = counts[::-1]

    # If nBases > 0, we have an easy case. Just add the bases.
    if nBases >= 0:
        if operations[0] == "M":
            counts[0] = counts[0] + nBases                      
        else:
            counts = [nBases] + counts
            operations = ["M"] + operations
    # If nBases < 0, we need to make sure we have room to delete the bases
    else:    
        for i in range(0,len(counts)):
           ct = counts[i]
           remainder = ct - abs(nBases)
           if remainder > 0:
                counts[i] = remainder
                break
           elif remainder == 0:
               counts[i] = ""
               operations[i] = ""
               break
           else:
               counts[i] = ""
               operations[i] = ""
               nBases = remainder

    if side == -1:
        operations = operations[::-1]
        counts = counts[::-1]

    result = ""
    for op, ct in zip(operations, counts):
        result = result + str(ct) + op
    return result
        
def dryRun_recordIndels(sam, outprefix, genome):
    """Records all insertions and deletions in the transcripts,
       but does not correct them """

    errorLog = outprefix + "_clean.TE.log"
    transcriptLog = outprefix + "_clean.log"
    eL = open(errorLog, 'w')
    tL = open(transcriptLog, 'w')

    eL.write("\t".join(["TranscriptID", "Position", "ErrorType", "Size", "Corrected", "ReasonNotCorrected"]) + "\n")
    tL.write("\t".join(["TranscriptID", "Mapping", \
                        "corrected_deletions", "uncorrected_deletions", "variant_deletions", \
                        "corrected_insertions", "uncorrected_insertions", "variant_insertions", \
                        "corrected_mismatches", "variant_mismatches", \
                        "corrected_NC_SJs", "uncorrected_NC_SJs"]) + "\n")
    spliceAnnot = {}
    with open(sam, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith("@"): # header line
                continue

            transcript = Transcript2(line, genome, spliceAnnot)
            
            if transcript.mapping == 0:
                logInfo = [transcript.QNAME, "unmapped", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
                tL.write("\t".join(logInfo) + "\n")
                continue
            if transcript.mapping == 2:
                logInfo = [transcript.QNAME, "non-primary", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]
                tL.write("\t".join(logInfo) + "\n")
                continue

            logInfo = [transcript.QNAME, "primary", "NA", 0, "NA", "NA", 0, "NA", "NA", "NA", "NA", "NA"]
            seqPos = 0
            genomePos = transcript.POS
            mergeOperations, mergeCounts = transcript.mergeMDwithCIGAR()

            for op,ct in zip(mergeOperations, mergeCounts):
                if op in ["M", "X"]:
                    seqPos += ct
                    genomePos += ct

                if op == "D":
                    ID = "_".join([transcript.CHROM, str(genomePos), str(genomePos + ct - 1)])
                    eL.write("\t".join([transcript.QNAME, ID, "Deletion", str(ct), "Uncorrected", "DryRun"]) + "\n")
                    logInfo[3] += 1
                    genomePos += ct
                   
                if op == "I":
                    ID = "_".join([transcript.CHROM, str(genomePos), str(genomePos + ct - 1)])
                    eL.write("\t".join([transcript.QNAME, ID, "Insertion", str(ct), "Uncorrected", "DryRun"]) + "\n")
                    logInfo[6] += 1
                    seqPos += ct
 
                if op == "S":
                    seqPos += ct

                if op in ["N", "H"]:
                    genomePos += ct
            tL.write("\t".join([str(i) for i in logInfo]) + "\n")
    eL.close()
    tL.close()
    return 

main()
