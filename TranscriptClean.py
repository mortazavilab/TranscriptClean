# TranscriptClean Version 2
# Dana Wyman, 11/20/2017
# In this version of TranscriptClean, mismatches and microindels in long reads are corrected in a SNP-aware fashion using the reference genome and a VCF file of whitelisted variants. Noncanonical splice junctions can also be corrected using a file of reference splice sites as long as the provided SAM file is splice aware.

from transcript2 import Transcript2
from spliceJunction import SpliceJunction
from intronBound import IntronBound
from optparse import OptionParser
import pybedtools
from pyfasta import Fasta
import os
#import re

def getOptions():
    parser = OptionParser()
    parser.add_option("--sam", dest = "sam", help = "Input file",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--genome", dest = "refGenome", help = "Reference genome fasta file. Should be the same one used during mapping to generate the sam file.",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--spliceJns", dest = "spliceAnnot", help = "Splice junction file obtained by mapping Illumina reads to the genome using STAR. More formats may be supported in the future. (Optional, but necessary for splice junction correction).", metavar = "FILE", type = "string", default = None)
    parser.add_option("--variants", dest = "variantFile",
                      help = "VCF formatted file of variants to avoid correcting away in the data (optional).", metavar = "FILE", type = "string", default = None )
    parser.add_option("--maxLenIndel", dest = "maxLenIndel",
                      help = "Maximum size indel to correct (Default: 5 bp)", type = "int", default = 5 )
    parser.add_option("--outprefix", dest = "outprefix",
                      help = "output file prefix. '_clean' plus a file extension will be added to the end.", metavar = "FILE", type = "string", default = "out")
    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    
    # Read in the reference genome. Treat coordinates as 0-based 
    print "Reading genome .............................."
    genome = Fasta(options.refGenome)

    if options.spliceAnnot != None:
        print "Processing annotated splice junctions ..."
        annotatedSpliceJns = processSpliceAnnotation(options.spliceAnnot)
    else:
        print "No splice annotation provided. Will skip splice junction correction."

    if options.variantFile != None:
        print "Processing variant file ................."
        #variants = processVariants(options.variantFile)
    else:
        print "No variant file provided. Transcript correction will not be SNP-aware."
        variants = {}

    print "Processing SAM file ........................."
    header, canTranscripts, noncanTranscripts = processSAM(options.sam, genome) 
    if len(noncanTranscripts) == 0: print "Note: No noncanonical transcripts found. This might mean that the sam file lacked the jM tag."

    print "Correcting mismatches and indels ............"
    correctMismatchesAndIndels(canTranscripts, genome, variants, options.maxLenIndel)
    correctMismatchesAndIndels(noncanTranscripts, genome, variants, options.maxLenIndel)

def processSAM(sam, genome):
    # This function extracts the SAM header (because we'll need that later) and creates a Transcript object for every sam transcript. 
    # Transcripts are returned two separate lists: one canonical and one noncanonical. 

    header = ""
    canTranscripts = {}
    noncanTranscripts = {}

    with open(sam, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                header = header + line + "\n"
                continue
            t = Transcript2(line, genome)

            # Filter out transcripts that are multimapping. The primary alignment will still be kept because its flag is <= 16
            if int(t.FLAG) > 16:
                continue

            # Filter out unmapped transcripts altogether
            if t.CHROM == "*":
                continue

            # Assign transcript to canonical or noncanonical group based on splice junctions
            if t.isCanonical == True:
                canTranscripts[t.QNAME] = t
            else:
                noncanTranscripts[t.QNAME] = t

    return header, canTranscripts, noncanTranscripts 

def processSpliceAnnotation(annotFile):
    # This function reads in the tab-separated STAR splice junction file and creates a bedtools object

    bedstr = ""
    o = open("TC-tmp.bed", 'w')
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
    o.close()

    # Convert bed file into BedTool object
    os.system('sort -k1,1 -k2,2n TC-tmp.bed > TC-tmp2.bed')
    bt = pybedtools.BedTool("TC-tmp2.bed")
    return bt

def correctMismatchesAndIndels(transcripts, genome, variants, maxLen):
    # This function corrects mismatches and indels up to size maxLen using the reference genome. If a variant file was provided, correction will be SNP-aware.
    
    for t in transcripts.keys():
        t = transcripts[t]

        origSeq = t.SEQ
        origCIGAR = t.CIGAR
        origNM = t.NM.split(":")[2]
        origMD = t.MD
        genomePos = t.POS

        print origCIGAR
        print origNM
        # Check for deletions, insertions, and mismatches. If none are present, we can skip this transcript
        if ("D" not in origCIGAR) and ("I" not in origCIGAR) and any(i in origNM for i in 'ACTGN') == False : continue
       
        newCIGAR = ""
        newSeq = ""
        MVal = 0
        seqPos = 0
        genomePos = t.POS
        
        cigarOperations, cigarCounts = t.splitCIGAR() 
        mdOperations, mdCounts = t.splitMD()        

        #for op,ct in zip(cigarOperations, cigarCounts):
        


        exit()

main()
