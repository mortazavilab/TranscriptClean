# TranscriptClean
# Author: Dana Wyman
# 3/1/2018
# -----------------------------------------------------------------------------
# Mismatches and microindels in long reads are corrected in a variant-aware 
# fashion using the reference genome and a VCF file of whitelisted variants. 
# Noncanonical splice junctions can also be corrected using a file of reference 
# splice sites.

from transcript2 import Transcript2
from spliceJunction import *
from intronBound import IntronBound
from optparse import OptionParser
import pybedtools
import dstruct
from pyfasta import Fasta
import os
import re
import copy

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

def prep_refs(orig_options):
    """ Process input files and store them in a reference dict """

    options = cleanup_options(orig_options)
    genomeFile = options.refGenome
    variantFile = options.variantFile
    sjFile = options.spliceAnnot
    outprefix = options.outprefix

    # Container for references
    refs = dstruct.Struct()

    # Read in the reference genome.
    print("Reading genome ..............................")
    refs["genome"] = Fasta(genomeFile)

    # Read in splice junctions
    if sjFile != None:
        print("Processing annotated splice junctions ...")
        refs.donors, refs.acceptors, refs.sjDict = processSpliceAnnotation(sjFile, outprefix)
    else:
        print("No splice annotation provided. Will skip splice junction correction.")
        refs.donors = None
        refs.acceptors = None
        refs["sjDict"] = {}

    # Read in variants
    if variantFile != None:
        print("Processing variant file .................")
        refs["snps"], refs["insertions"], refs["deletions"] = processVCF(variantFile, options.maxLenIndel)
    else:
        print("No variant file provided. Transcript correction will not be variant-aware.")
        refs["snps"] = refs["insertions"] = refs["deletions"] = {}

    return options, refs

def cleanup_options(options):
    """ Clean up input options by casting to appropriate types etc. """
    options.maxLenIndel = int(options.maxLenIndel)
    options.maxSJOffset = int(options.maxSJOffset)
    options.indelCorrection = (options.correctIndels).lower()
    options.mismatchCorrection = (options.correctMismatches).lower()
    options.sjCorrection = (options.correctSJs).lower()

    return options

def setup_outfiles(options, process = ""):
    """ Set up output files. If running in parallel, label with a process ID """

    outfiles = dstruct.Struct()
    outprefix = options.outprefix
    
    # Open sam, fasta, and log  outfiles
    oSam = open(outprefix + "_clean" + process + ".sam", 'w')
    oFa = open(outprefix + "_clean" + process + ".fa", 'w')
    transcriptLog = open(outprefix + "_clean" + process + ".log", 'w')
    transcriptErrorLog = open(options.outprefix + "_clean" + process + ".TE.log", 'w')

    # Add headers to logs
    transcriptLog.write("\t".join(["TranscriptID", "Mapping", 
                        "corrected_deletions", "uncorrected_deletions", 
                        "variant_deletions", "corrected_insertions", 
                        "uncorrected_insertions", "variant_insertions", \
                        "corrected_mismatches", "uncorrected_mismatches", \
                        "corrected_NC_SJs", "uncorrected_NC_SJs"]) + "\n")

    transcriptErrorLog.write("\t".join(["TranscriptID", "Position", "ErrorType", 
                             "Size", "Corrected", "ReasonNotCorrected"]) + "\n")
    
    outfiles.sam = oSam
    outfiles.fasta = oFa
    outfiles.log = transcriptLog
    outfiles.TElog = transcriptErrorLog

    return outfiles

def close_outfiles(outfiles):
    """ Close all of the outout files """

    for f in list(outfiles.values()):
        f.close()

    return

def transcript_init(transcript_line, options, refs, outfiles):
    """ Attempt to initialize a Transcript object from a SAM entry. If the
        transcript alignment is unmapped or non-primary, then output that 
        information to the outfiles, but return None. Otherwise, return the
        transcript object. """

    outSam = outfiles.sam
    outFa = outfiles.fasta

    if transcript_line.startswith("@"): # header line
        outSam.write(transcript_line + "\n")
        return None
    
    # Init transcript object and log entry
    transcript = Transcript2(transcript_line, refs.genome, refs.sjDict)
    logInfo = init_log_info()
    logInfo.TranscriptID = transcript.QNAME

    # Unmapped/multimapper cases
    if transcript.mapping != 1:
        if transcript.mapping == 0:
            logInfo.Mapping = "unmapped"
        elif transcript.mapping == 2:
            logInfo.Mapping = "non-primary"

        # Only output the transcript if Primary Only option is off
        if primaryOnly == "false":
            outSam.write(transcript_line + "\n")
            outFa.write(Transcript2.printableFa(transcript) + "\n")

        return None, logInfo
    else:
        logInfo.Mapping = "primary"
        return transcript, logInfo 


def correct_transcript(transcript_line, options, refs, outfiles):
    """ Given a line from a SAM file, create a transcript object. If it's a
        primary alignment, then perform the corrections specified in the 
        options.
    """
    outSam = outfiles.sam
    outFa = outfiles.fasta

    if transcript_line.startswith("@"): # header line
        outSam.write(transcript_line + "\n")
        return

    transcript, logInfo = transcript_init(transcript_line, options, refs, outfiles)
    
    # Correct the transcript 
    if transcript != None:
        print(transcript.QNAME)

        # Mismatch correction
        if options.mismatchCorrection == "true":
            correctMismatches(transcript, refs.genome, refs.snps, 
                              logInfo, outfiles.TElog)
        
        if options.indelCorrection == "true":
            # Insertion correction
            correctInsertions(transcript, refs.genome, refs.insertions,
                              options.maxLenIndel, logInfo, outfiles.TElog)

            # Deletion correction
            correctDeletions(transcript, refs.genome, refs.deletions,
                              options.maxLenIndel, logInfo, outfiles.TElog)

        # NCSJ correction
        if refs.sjDict != {} and transcript.isCanonical == False:
            transcript = cleanNoncanonical(transcript, refs, options.maxSJOffset, 
                                           logInfo, outfiles.TElog)
    
    # Output transcript log entry 
    write_to_transcript_log(logInfo, outfiles.log)

    # Write transcript to sam and fasta file
    outSam.write(transcript.printableSAM() + "\n")
    outFa.write(transcript.printableFa() + "\n")
    return    

def validate_chroms(genome, sam):
    """ Make sure that every chromosome in the SAM file also exists in the
        reference genome. This is a common source of crashes """

    fasta_chroms = set(genome.keys())
    sam_chroms = set()
    with open(sam, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                continue
            else:
                chrom = line.split("\t")[2]
                sam_chroms.add(chrom)

    # Check whether all of the sam chromosomes are in the fasta file. 
    # If not, raise an error
    if not sam_chroms.issubset(fasta_chroms):
        sam_chroms = "{" + ", ".join(['"' + str(x) + '"' for x in sam_chroms]) + '}'
        fasta_chroms = "{" + ", ".join(['"' + str(x) + '"' for x in fasta_chroms]) + '}'
        print(sam_chroms)
        print(fasta_chroms)
        error_msg = "One or more SAM chromosomes were not found in the " +\
                    "fasta reference.\n" + \
                    "SAM chromosomes:\n" + sam_chroms + "\n" + \
                    "FASTA chromosomes:\n" + fasta_chroms + "\n" + \
                    "One common cause is when the fasta headers contain more than " +\
                    "one word. If this is the case, try trimming the headers " + \
                    "to only the chromosome name (i.e. '>chr1')."
        raise ValueError(error_msg)
    return

def main():
    orig_options = getOptions()
    options, refs = prep_refs(orig_options)
    validate_chroms(refs.genome, options.sam)
    outfiles = setup_outfiles(options)
    samFile = options.sam
    
    if options.dryRun == True:
        dryRun(samFile, outfiles, refs)
        close_outfiles(outfiles)
        return

    print("Processing transcripts in SAM file .........................")

    with open(samFile, 'r') as f:
        for transcript_line in f:
            correct_transcript(transcript_line, options, refs, outfiles)            

    close_outfiles(outfiles)


def processSpliceAnnotation(annotFile, outprefix):
    """ Reads in the tab-separated STAR splice junction file and creates a 
        bedtools object. Also creates a dict (annot) to allow easy lookup 
        to find out if a splice junction is annotated or not """

    bedstr = ""
    annot = set()
    donor_file = outprefix + "_ref_splice_donors_tmp.bed"
    acceptor_file = outprefix + "_ref_splice_acceptors_tmp.bed"
    o_donor = open(donor_file, 'w')
    o_acceptor = open(acceptor_file, 'w')
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

            # ID the splice donor/acceptor identity
            if strand == "+":
                file1 = o_donor
                file2 = o_acceptor
                type1 = "donor"
                type2 = "acceptor"
            elif strand == "-":
                file1 = o_acceptor
                file2 = o_donor
                type1 = "acceptor"
                type2 = "donor"

            # Make one bed entry for each end of the junction and write to
            # splice donor and acceptor files
            bed1 = "\t".join([chrom, str(start - 1), str(start), ".", uniqueReads, strand])
            bed2 = "\t".join([chrom, str(end - 1), str(end), ".", uniqueReads, strand])
            file1.write(bed1 + "\n")
            file2.write(bed2 + "\n")

            annot.add("_".join([chrom, str(start), strand, type1]))
            annot.add("_".join([chrom, str(end), strand, type2]))
    o_donor.close()
    o_acceptor.close()

    # Convert bed files into BedTool objects
    donor_sorted = outprefix + "_ref_splice_donors_tmp.sorted.bed"
    acceptor_sorted = outprefix + "_ref_splice_acceptors_tmp.sorted.bed"
    os.system('bedtools sort -i ' + donor_file + ' >' + donor_sorted)
    os.system('bedtools sort -i ' + acceptor_file + ' >' + acceptor_sorted)
    #spliceJnBedTool = pybedtools.BedTool(fNameSorted)
    splice_donor_bedtool = pybedtools.BedTool(donor_sorted)
    splice_acceptor_bedtool = pybedtools.BedTool(acceptor_sorted)
    os.system("rm " + donor_file)
    os.system("rm " + acceptor_file)
    return splice_donor_bedtool, splice_acceptor_bedtool, annot

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

def correctInsertions(transcript, genome, variants, maxLen, logInfo, eL):
    """ Corrects insertions up to size maxLen using the reference genome. 
        If a variant file was provided, correction will be SNP-aware. """

    logInfo.uncorrected_insertions = 0
    logInfo.corrected_insertions = 0

    origSeq = transcript.SEQ
    origCIGAR = transcript.CIGAR
    transcript_ID = transcript.QNAME    
 
    cigarOps,cigarCounts = transcript.splitCIGAR() 

    # Check for insertions. If none are present, we can skip this transcript
    if "I" not in origCIGAR: return

    newCIGAR = ""
    newSeq = ""
    MVal = 0
    seqPos = 0

    # Start at position in the genome where the transcript starts.
    genomePos = transcript.POS 

    # Iterate over operations to sequence and repair insertions
    for op,ct in zip(cigarOps, cigarCounts):

        currPos = transcript.CHROM + ":" + str(genomePos) + "-" + \
                  str(genomePos + ct - 1)          

        if op == "M":
             newSeq = newSeq + origSeq[seqPos:seqPos + ct]
             MVal += ct
             seqPos += ct
             genomePos += ct
         
        if op == "I":
            ID = "_".join([transcript.CHROM, str(genomePos), 
                           str(genomePos + ct - 1)])

            # Only insertions of a given size are corrected
            if ct <= maxLen:
                # Check if the insertion is in the optional variant catalog.
                if ID in variants:
                    # The insertion perfectly matches a variant position. 
                    # Leave the sequence alone if it matches an allele sequence.
                    currSeq = origSeq[seqPos:seqPos + ct]
                    if currSeq in variants[ID]:
                        logInfo.variant_insertions += 1
                        #Transcript2.addVariantInsertion(t)
                        errorEntry = "\t".join([transcript_ID, ID, "Insertion", 
                                                str(ct), "Uncorrected", 
                                                "VariantMatch"])
                        eL.write(errorEntry + "\n")  

                        # Leave insertion in 
                        MVal, newCIGAR = endMatch(MVal, newCIGAR)
                        newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                        newCIGAR = newCIGAR + str(ct) + op
                        seqPos += ct
                        continue

                # Correct insertion
                errorEntry = "\t".join([transcript_ID, ID, "Insertion", str(ct), 
                                        "Corrected", "NA"])
                logInfo.corrected_insertions += 1
                eL.write(errorEntry + "\n")

                # Subtract the inserted bases by skipping them. 
                # GenomePos stays the same, as does MVal
                seqPos += ct
            else: # Move on without correcting insertion because it is too big
                errorEntry = "\t".join([transcript_ID, ID, "Insertion", str(ct), 
                                        "Uncorrected", "TooLarge"])
                logInfo.uncorrected_insertions += 1
                eL.write(errorEntry + "\n")
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
    transcript.CIGAR = newCIGAR
    transcript.SEQ = newSeq
        
    return    

def correctDeletions(transcript, genome, variants, maxLen, logInfo, eL):
    """ Corrects deletions up to size maxLen using the reference genome. 
        If a variant file was provided, correction will be variant-aware."""

    logInfo.uncorrected_deletions = 0
    logInfo.corrected_deletions = 0

    transcript_ID = transcript.QNAME
    chrom = transcript.CHROM
    origSeq = transcript.SEQ
    origCIGAR = transcript.CIGAR

    cigarOps,cigarCounts = transcript.splitCIGAR()

    # Check for deletions. If none are present, we can skip this transcript
    if "D" not in origCIGAR: 
        return

    newCIGAR = ""
    newSeq = ""
    MVal = 0
    seqPos = 0

    # Start at position in the genome where the transcript starts.
    genomePos = transcript.POS

    # Iterate over operations to sequence and repair mismatches and microindels
    for op,ct in zip(cigarOps, cigarCounts):

        currPos = chrom + ":" + str(genomePos) + "-" + \
                  str(genomePos + ct - 1) 

        if op == "M":
             newSeq = newSeq + origSeq[seqPos:seqPos + ct]
             MVal += ct
             seqPos += ct
             genomePos += ct

        if op == "D":
            ID = "_".join([chrom, str(genomePos), str(genomePos + ct - 1)])
            if ct <= maxLen:

                # Check if the deletion is in the optional variant catalog.
                if ID in variants:
                    # The deletion perfectly matches a deletion variant. Leave the deletion in.
                    currSeq = genome.sequence({'chr': chrom, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                    errorEntry = "\t".join([transcript_ID, ID, "Deletion", str(ct), "Uncorrected", "VariantMatch"])
                    logInfo.variant_deletions += 1
                    eL.write(errorEntry + "\n")

                    MVal, newCIGAR = endMatch(MVal, newCIGAR)
                    genomePos += ct
                    newCIGAR = newCIGAR + str(ct) + op
                    continue

                # Correct deletion if we're not in variant-aware mode
                errorEntry = "\t".join([transcript_ID, ID, "Deletion", str(ct), "Corrected", "NA"])
                logInfo.corrected_deletions += 1
                eL.write(errorEntry + "\n")

                # Add the missing reference bases
                refBases = genome.sequence({'chr': chrom, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                newSeq = newSeq + refBases
                genomePos += ct
                MVal += ct

            # Deletion is too big to fix
            else:
                errorEntry = "\t".join([transcript_ID, ID, "Deletion", str(ct), "Uncorrected", "TooLarge"])
                logInfo.uncorrected_deletions += 1
                eL.write(errorEntry + "\n")

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

    transcript.CIGAR = newCIGAR
    transcript.SEQ = newSeq
    transcript.NM, transcript.MD = transcript.getNMandMDFlags(genome)
    return 


def correctMismatches(transcript, genome, variants, logInfo, eL):
    """ This function corrects mismatches in the provided transcript. If a 
        variant file was provided, correction will be variant-aware."""
   
    logInfo.uncorrected_mismatches = 0
    logInfo.corrected_mismatches = 0

    origSeq = transcript.SEQ
    origCIGAR = transcript.CIGAR
    origMD = transcript.MD

    # Check for mismatches. If none are present, we can skip this transcript
    if any(i in origMD.upper() for i in 'ACTGN') == False : 
        return
        
    newCIGAR = ""
    newSeq = ""
    MVal = 0
    seqPos = 0
    genomePos = transcript.POS

    # Merge CIGAR and MD tag information so that we have the locations of all 
    # insertions, deletions, and mismatches
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
                    errorEntry = "\t".join([transcript.QNAME, ID, "Mismatch", 
                                            str(ct), "Uncorrected", 
                                            "VariantMatch"])
                    eL.write(errorEntry + "\n")
                    logInfo.uncorrected_mismatches += 1

                    # Keep the base as-is
                    newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                    MVal += ct
                    seqPos += ct
                    genomePos += ct
                    continue
            # Otherwise, correct the mismatch to reference base
            errorEntry = "\t".join([transcript.QNAME, ID, "Mismatch", 
                                    str(ct), "Corrected", "NA"])
            eL.write(errorEntry + "\n")
            logInfo.corrected_mismatches += 1

            # Change sequence base to the reference base at this position
            newSeq = newSeq + genome.sequence({'chr': transcript.CHROM, 
                                               'start': genomePos, 
                                               'stop': genomePos + ct - 1}, 
                                               one_based=True)
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

def cleanNoncanonical(transcript, refs, maxDist, logInfo, TElog):
   
    logInfo.corrected_NC_SJs = 0
    logInfo.uncorrected_NC_SJs = 0 

    if transcript.isCanonical == True:
        return transcript

    # Iterate over all junctions and attempt to correct
    for splice_jn_num in range(0,len(transcript.spliceJunctions)):
        junction = transcript.spliceJunctions[splice_jn_num]
        if junction.isCanonical:
            continue
        else:
            temp_transcript = copy.copy(transcript)
            ID = "_".join([junction.chrom, str(junction.start), str(junction.end)]) 
            correction_status, reason, dist = attempt_jn_correction(temp_transcript, 
                                                                    splice_jn_num,
                                                                    refs.genome,
                                                                    refs.donors,
                                                                    refs.acceptors,
                                                                    refs.sjDict,
                                                                    maxDist)

            # If there were no problems during correction, replace the original 
            # transcript with the modified copy. 
            if correction_status == True:
                transcript = temp_transcript
                logInfo.corrected_NC_SJs += 1
                status = "Corrected"
            else:
                logInfo.uncorrected_NC_SJs += 1
                status = "Uncorrected"
                
            # Update transcript error log
            errorEntry = "\t".join([transcript.QNAME, ID, "NC_SJ_boundary",
                                    str(dist), status, reason])
            TElog.write(errorEntry + "\n")

    return transcript

def find_closest_bound(sj_bound, ref_bounds):
    """ Given one side of a splice junction, find the closest reference """

    # Create a Bedtool object for the bound
    bed_pos = pybedtools.BedTool(sj_bound.getBED(), from_string=True)

    # Run Bedtools Closest operation
    closest = bed_pos.closest(ref_bounds, s=True, D="ref", t="first")[0]

    # Create an object to represent the closest match
    # Coordinates are 0-based since they are coming from BED
    obj_closest = dstruct.Struct()
    obj_closest.chrom = closest[6]
    obj_closest.start = int(closest[7])
    obj_closest.end = int(closest[8])
    obj_closest.dist = int(closest[-1])
     
    return obj_closest

def find_closest_ref_junction(junction, ref_donors, ref_acceptors):
    """ Given a splice junction object and a splice reference, locate the 
        closest junction in the provided splice reference BedTool object"""

    # Get the splice donor and acceptor of the junction
    donor = junction.get_splice_donor()
    acceptor = junction.get_splice_acceptor()

    # Find the closest match for each
    closest_ref_donor = find_closest_bound(donor, ref_donors)
    closest_ref_acceptor = find_closest_bound(acceptor, ref_acceptors)

    return closest_ref_donor, closest_ref_acceptor


def update_post_ncsj_correction(transcript, splice_jn_num, genome, spliceDict):
    """ After correcting a noncanonical splice junction, perform the 
        following updates:
        - Reassign the position of the junction (based on its bounds)
        - Recompute the splice motif and assign to junction
        - Recompute NM/MD tags
    """
    junction = transcript.spliceJunctions[splice_jn_num]
    junction.recheckPosition()
    junction.checkSpliceMotif(genome, spliceDict)
    transcript.NM, transcript.MD = transcript.getNMandMDFlags(genome)
    transcript.jM, transcript.jI = transcript.get_jM_jI_tags_from_sjs() 
    transcript.isCanonical = transcript.recheckCanonical()
    return

def attempt_jn_correction(transcript, splice_jn_num, genome, ref_donors, 
                          ref_acceptors, spliceDict, maxDist):
    """ Given a noncanonical splice junction, try to correct it by locating 
        a nearby annotated junction within the allowable distance.

        Returns:
            Corrected (True/False)
            Reason for no correction ("NA"/"TooFarFromAnnotJn"/"Other")
            CombinedDist
     """

    junction = transcript.spliceJunctions[splice_jn_num]
    
    # Find the closest reference junction
    ref_donor, ref_acceptor = find_closest_ref_junction(junction, ref_donors, 
                                                         ref_acceptors)

    # Compute the overall distance of the original junction from the reference
    combined_dist = combinedJunctionDist(ref_donor.dist, ref_acceptor.dist)

    # Only attempt to rescue junction boundaries that are within maxDist bp of 
    # an annotated junction
    if combined_dist > maxDist or abs(ref_donor.dist) > 2*maxDist or \
                                  abs(ref_acceptor.dist) > 2*maxDist:
        return False, "TooFarFromAnnotJn", combined_dist
   
    try:
        # Attempt to fix the splice donor side
        donor = junction.get_splice_donor()
        transcript.SEQ, transcript.CIGAR = fix_one_side_of_junction(transcript.CHROM, 
                                           transcript.POS, splice_jn_num, 
                                           donor, ref_donor.dist, genome, 
                                           transcript.SEQ, transcript.CIGAR)

        # Attempt to fix the splice acceptor side
        acceptor = junction.get_splice_acceptor()
        transcript.SEQ, transcript.CIGAR = fix_one_side_of_junction(transcript.CHROM,
                                           transcript.POS, splice_jn_num,
                                           acceptor, ref_acceptor.dist, genome,
                                           transcript.SEQ, transcript.CIGAR)
        # Now, perform updates:
        update_post_ncsj_correction(transcript, splice_jn_num, genome, spliceDict)

    except:
        return False, "Other", combined_dist

    return True, "NA", combined_dist 
    

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


def group_CIGAR_by_exon_and_intron(CIGAR, seq):
    """ Given a CIGAR string, group all of the operations by exon or intron, 
        resulting in one string per exon/intron. Also, subdivide the provided
        sequence by exon as well."""
    
    # Inits
    currExonStr = ""
    currExonCIGAR = ""
    currSeqIndex = 0
    operations, counts = splitCIGAR(CIGAR)

    # Outputs
    exonSeqs = []
    exonCIGARs = []
    intronCIGARs = []

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
    # Append last exon entries
    exonSeqs.append(currExonStr)
    exonCIGARs.append(currExonCIGAR)

    return exonCIGARs, intronCIGARs, exonSeqs

def fix_one_side_of_junction(chrom, transcript_start, jn_number, intronBound, d, 
                             genome, seq, oldCIGAR):
    """ Corrects one side of a noncanonical splice junction by the specified 
        number of bases. This involves modifying the sequence and CIGAR string.

        Returns modified CIGAR string and sequence
    """

    # If d = 0, then the bound already matches the annotation. 
    if d == 0:
        return seq, oldCIGAR

    # First, split CIGAR by exon/intron, as well as the sequence
    exonCIGARs, intronCIGARs, exonSeqs = group_CIGAR_by_exon_and_intron(oldCIGAR, seq)

    # Now go and fix the specific splice junction piece
    if intronBound.bound == 0:
        targetExon = jn_number
        exon = exonSeqs[targetExon]
        if d > 0: # Need to add d bases from reference to end of exon. Case 1.
            # For CIGAR string,
            exonEnd = intronBound.pos - 1
            seqIndex = exonEnd - transcript_start + 1
            refAdd = genome.sequence({'chr': chrom, 'start': exonEnd + 1, 'stop': exonEnd + d}, one_based=True)
            exonSeqs[targetExon] = exon + refAdd
            intronBound.pos += d
        if d < 0: # Need to subtract from end of exon sequence. Case 3
            exonSeqs[targetExon] = exon[0:d]
            intronBound.pos += d
        intronCIGARs[jn_number] -= d
        exonCIGARs[targetExon] = editExonCIGAR(exonCIGARs[targetExon], -1, d)
    else:
        targetExon = jn_number + 1
        exon = exonSeqs[targetExon]
        if d < 0: # Need to add d bases from reference to start of exon sequence. Case 2.
            exonStart = intronBound.pos + 1
            seqIndex = exonStart - transcript_start + 1
            refAdd = genome.sequence({'chr': chrom, 'start': exonStart - abs(d), 'stop': exonStart - 1}, one_based=True)
            exonSeqs[targetExon] = refAdd + exon
            intronBound.pos += d
        if d > 0: # Need to subtract from start of exon sequence. Case 4
            exonSeqs[targetExon] = exon[d:]
            intronBound.pos += d

        # Modify exon string
        exonCIGARs[targetExon] = editExonCIGAR(exonCIGARs[targetExon], 0, -d)
        intronCIGARs[jn_number] += d

    newSeq = str(''.join(exonSeqs))

    # Paste together the new CIGAR string
    newCIGAR = ""
    for i in range(0,len(intronCIGARs)):
        newCIGAR = newCIGAR + exonCIGARs[i] + str(intronCIGARs[i]) + "N"
    newCIGAR = newCIGAR + exonCIGARs[-1]

    return newSeq, newCIGAR


def check_CIGAR_validity(CIGAR):
    """ Returns 'True' if CIGAR string passes tests, and "False" otherwise.
        Checks to make sure that 
        1) Introns (N) are only ever followed by M operation
        2) Introns never follow an insertion/deletion
    """

    if "-" in CIGAR:
        return False

    operations, counts = splitCIGAR(CIGAR)
    prev = ""
    for op, ct in zip(operations, counts):
        if prev == "N" and op != "M":
            return False
        if prev == "D" and op == "N":
            return False
        
        prev = op
    return True
 
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
        
def init_log_info():
    """ Initialize a log_info struct to track the error types and correction
        status for a single transcript """

    logInfo = dstruct.Struct()
    logInfo.TranscriptID = None
    logInfo.Mapping = None
    logInfo.corrected_deletions = "NA"
    logInfo.uncorrected_deletions = "NA"
    logInfo.variant_deletions = "NA"
    logInfo.corrected_insertions = "NA"
    logInfo.uncorrected_insertions = "NA"
    logInfo.variant_insertions = "NA"
    logInfo.corrected_mismatches = "NA"
    logInfo.uncorrected_mismatches = "NA"
    logInfo.corrected_NC_SJs = "NA"
    logInfo.uncorrected_NC_SJs = "NA"

    return logInfo

def write_to_transcript_log(logInfo, tL):
    """ Write a transcript log entry to output """

    log_strings = [ str(x) for x in [logInfo.TranscriptID, logInfo.Mapping,
                                     logInfo.corrected_deletions,
                                     logInfo.uncorrected_deletions,
                                     logInfo.variant_deletions,
                                     logInfo.corrected_insertions,
                                     logInfo.uncorrected_insertions,
                                     logInfo.variant_insertions,
                                     logInfo.corrected_mismatches,
                                     logInfo.uncorrected_mismatches,
                                     logInfo.corrected_NC_SJs,
                                     logInfo.uncorrected_NC_SJs] ]

    tL.write("\t".join(log_strings) + "\n")
    return


def dryRun(sam, outfiles, refs):
    """Records all mismatches, insertions, and deletions in the transcripts,
       but does not correct them """

    print("Dry run mode: Cataloguing indels and mismatches.........")
    tL = outfiles.log
    eL = outfiles.TElog
    genome = refs.genome
    spliceAnnot = {}

    with open(sam, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith("@"): # header line
                continue

            # Set up log entry for transcript
            transcript = Transcript2(line, genome, spliceAnnot)
            logInfo = init_log_info()
            logInfo.TranscriptID = transcript.QNAME
            
            if transcript.mapping == 0:
                logInfo.Mapping = "unmapped"
                write_to_transcript_log(logInfo, tL)
                continue
            if transcript.mapping == 2:
                logInfo.Mapping = "non-primary"
                write_to_transcript_log(logInfo, tL)
                continue

            # For primary alignments, modify logInfo
            logInfo.Mapping = "primary" 
            logInfo.uncorrected_deletions = 0
            logInfo.uncorrected_insertions = 0
            logInfo.uncorrected_mismatches = 0

            # Iterate over CIGAR to catalogue indels and mismatches
            seqPos = 0
            genomePos = transcript.POS
            mergeOperations, mergeCounts = transcript.mergeMDwithCIGAR()

            for op,ct in zip(mergeOperations, mergeCounts):
                if op == "M":
                    seqPos += ct
                    genomePos += ct
               
                if op == "X":
                    ID = "_".join([transcript.CHROM, str(genomePos), 
                                   str(genomePos + ct - 1)])
                    eL.write("\t".join([transcript.QNAME, ID, "Mismatch", 
                                        str(ct), "Uncorrected", "DryRun"]) + "\n")
                    logInfo.uncorrected_mismatches += 1
                    seqPos += ct
                    genomePos += ct

                if op == "D":
                    ID = "_".join([transcript.CHROM, str(genomePos), 
                                   str(genomePos + ct - 1)])
                    eL.write("\t".join([transcript.QNAME, ID, "Deletion", 
                                        str(ct), "Uncorrected", "DryRun"]) + "\n")
                    logInfo.uncorrected_deletions += 1
                    genomePos += ct
                   
                if op == "I":
                    ID = "_".join([transcript.CHROM, str(genomePos), 
                                   str(genomePos + ct - 1)])
                    eL.write("\t".join([transcript.QNAME, ID, "Insertion", 
                                        str(ct), "Uncorrected", "DryRun"]) + "\n")
                    logInfo.uncorrected_insertions += 1
                    seqPos += ct
 
                if op == "S":
                    seqPos += ct

                if op in ["N", "H"]:
                    genomePos += ct
            write_to_transcript_log(logInfo, tL)
    return 

if __name__ == '__main__':
    main()
