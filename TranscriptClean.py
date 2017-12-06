# TranscriptClean Version 2
# Dana Wyman, 11/20/2017
# In this version of TranscriptClean, mismatches and microindels in long reads are corrected in a SNP-aware fashion using the reference genome and a VCF file of whitelisted variants. Noncanonical splice junctions can also be corrected using a file of reference splice sites as long as the provided SAM file is splice aware.

## TODO: Add a line to the sam header that explains which input options TranscriptClean was run with
## TODO: Add data structure/system to track transcript changes
## TODO: Add options input arguments to run only certain parts of TranscriptClean

from transcript2 import Transcript2
from spliceJunction import SpliceJunction
from intronBound import IntronBound
from optparse import OptionParser
import pybedtools
from pyfasta import Fasta
import os
import re

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
    parser.add_option("--maxSJOffset", dest = "maxSJOffset",
                      help = "Maximum distance from annotated splice junction to correct (Default: 5 bp)", type = "int", default = 5 )
    parser.add_option("--outprefix", dest = "outprefix",
                      help = "output file prefix. '_clean' plus a file extension will be added to the end.", metavar = "FILE", type = "string", default = "out")
    # Options that control which sections to run
    parser.add_option("--correctMismatches", dest = "correctMismatches",
                      help = "If set to false, TranscriptClean will skip mismatch correction. Default: True", type = "string", default = "true" )
    parser.add_option("--correctIndels", dest = "correctIndels",
                      help = "If set to false, TranscriptClean will skip indel correction. Default: True", type = "string", default = "true")
    

    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()

    # Read in the reference genome. 
    print "Reading genome .............................."
    genome = Fasta(options.refGenome)

    if options.spliceAnnot != None:
        print "Processing annotated splice junctions ..."
        annotatedSpliceJns = processSpliceAnnotation(options.spliceAnnot)
    else:
        print "No splice annotation provided. Will skip splice junction correction."

    if options.variantFile != None:
        print "Processing variant file ................."
        snps, insertions, deletions = processVCF(options.variantFile, options.maxLenIndel)
    else:
        print "No variant file provided. Transcript correction will not be SNP-aware."
        snps = {}
        insertions = {}
        deletions = {}

    print "Processing SAM file ........................."
    header, canTranscripts, noncanTranscripts = processSAM(options.sam, genome) 
    if len(noncanTranscripts) == 0: print "Note: No noncanonical transcripts found. If this is unexpected, check whether the sam file lacks the jM tag."
     
    if options.correctMismatches.lower() == "true":
        print "Correcting mismatches (canonical transcripts)............"
        correctMismatches(canTranscripts, genome, snps)
     
    if options.correctIndels.lower() == "true":
        print "Correcting insertions (canonical transcripts)............"
        correctInsertions(canTranscripts, genome, insertions, options.maxLenIndel)
        print "Correcting deletions (canonical transcripts)............"
        correctDeletions(canTranscripts, genome, deletions, options.maxLenIndel)

    if len(noncanTranscripts) > 0:
        if options.correctMismatches.lower() == "true":
            print "Correcting mismatches (noncanonical transcripts)............"
            correctMismatches(noncanTranscripts, genome, snps)
        if options.correctIndels.lower() == "true":
            print "Correcting insertions (noncanonical transcripts)............"
            correctInsertions(noncanTranscripts, genome, snps, options.maxLenIndel)
            print "Correcting deletions (noncanonical transcripts)............"
            correctDeletions(noncanTranscripts, genome, snps, options.maxLenIndel)

        print "Rescuing noncanonical junctions............."
        cleanNoncanonical(noncanTranscripts, annotatedSpliceJns, genome, options.maxSJOffset)

    print "Writing output to sam file and fasta file.................."

    # Generate the output files
    oSam = open(options.outprefix + "_clean.sam", 'w')
    oFa = open(options.outprefix + "_clean.fa", 'w')
    oLog = open(options.outprefix + "_clean.log", 'w')
    oSam.write(header)
    writeTranscriptOutput(canTranscripts, oSam, oFa, oLog, genome)
    writeTranscriptOutput(noncanTranscripts, oSam, oFa, oLog, genome)

    oSam.close()
    oFa.close()

def writeTranscriptOutput(transcripts, outSam, outFa, outLog, genome):

    for t in transcripts.keys():
        print t
        currTranscript = transcripts[t]
        outSam.write(Transcript2.printableSAM(currTranscript, genome) + "\n")
        outFa.write(Transcript2.printableFa(currTranscript) + "\n")
        outLog.write(t + "\t" + ";".join(currTranscript.log) + "\n")
    return

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

def processVCF(vcf, maxLen):
    # This function reads in variants from a VCF file and stores them.
    # SNPs are simple- they get stored in a dictionary by position.
    # Insertions get placed in a temporary bed file, as do deletions.
    SNPs = {}
    insertions = {}
    deletions = {}
    #insertionFile = open("tmp_vcf_insertions.bed", 'w')
    #deletionFile = open("tmp_vcf_deletions.bed", 'w')

    with open(vcf, 'r') as f:
        for line in f:
            line = line.strip()
    
            if line.startswith("#"): continue
            fields = line.split()
            
            chrom = "chr" + fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(",")
 
            refLen = len(ref)

            # It is possible for a single position to have more than one type of alternate allele (ie mismatch or indel). So it is necessary to treat each alternate allele separately.
            for allele in alt:
                altLen = len(allele)
                
                # SNP case
                if refLen == altLen == 1:
                    ID = chrom + "_" + str(pos)
                    if ID not in SNPs:
                        SNPs[ID] = [allele]
                    else:
                        SNPs[ID] = SNPs[ID].append(allele)
                # Deletion
                if refLen > altLen:
                    delSize = refLen - altLen
                    if delSize > maxLen: continue
                    # VCF files are one-based. Inserton/deletion sequences always start with the preceding normal reference base, so we need to add one to pos
                    # This principle holds for both insertions and deletions.
                    pos = pos + 1
                    ID = "_".join([ chrom, str(pos), str(pos + delSize - 1)])
                    deletions[ID] = [allele] 
                # Insertion
                if refLen < altLen:
                    inSize = altLen - refLen
                    if inSize > maxLen: continue
                    pos = pos + 1
                    ID = "_".join([ chrom, str(pos), str(pos + inSize - 1)])
                    insertions[ID] = [allele]
                # Deletion
                #if refLen > altLen:
                #    delSize = refLen - altLen
                    # VCF files are one-based. However, we don't need to subtract 1 from the start position when converting to 0-based BED here
                    # for the simple reason that in VCFs, inserton/deletion sequences always start with the preceding normal reference base.
                    # This principle holds for both insertions and deletions. 
                #    deletionFile.write("\t".join([ chrom, str(pos), str(pos + delSize), ".", "+", "."]) + "\n")
                    
                # Insertion
                #if refLen < altLen:
                #    inSize = altLen - refLen
                #    insertionFile.write("\t".join([ chrom, str(pos), str(pos + inSize), ".", "+", "."]) + "\n")

    #insertionFile.close()
    #deletionFile.close()

    #os.system('sort -k1,1 -k2,2n tmp_vcf_insertions.bed > sorted_tmp_vcf_insertions.bed')
    #os.system('sort -k1,1 -k2,2n tmp_vcf_deletions.bed > sorted_tmp_vcf_deletions.bed')
    #insertions = pybedtools.BedTool("sorted_tmp_vcf_insertions.bed")
    #deletions = pybedtools.BedTool("sorted_tmp_vcf_deletions.bed")
    return SNPs, insertions, deletions

def intersectWithVariants(transcriptIndels, variants):
    # This function intersects all transcript insertions or deletions with the variant set to look for the closest overlap. It returns a dictionary where each key chrom_start_end points to the extent of overlap.

    result = {}
    intersection = str(transcriptIndels.intersect(variants, wao = True)).split("\n")
    for line in intersection:
        info = line.split()
        if len(info) == 0: continue
        if info[-1] == "0": continue
        
        print line
        chrom = info[0]
        start = str(int(info[1]) + 1) #(convert back to 1-based)
        end = info[2]
        trID = info[3]
        overlap = int(info[-1])
        ID = "_".join([chrom, start, end, trID])

        # Add the ID to the result, along with the amount of overlap. Only keep the match with the largest overlap if there is more than one for a transcript.
        if ID not in result:
            result[ID] = overlap
        else:
            if overlap > result[ID]:
                result[ID] = overlap

    return result                

def correctInsertions(transcripts, genome, variants, maxLen):
    # This function corrects insertions up to size maxLen using the reference genome. If a variant file was provided, correction will be SNP-aware.

    # If variants are provided, get all transcript insertions and intersect them with the variant catalog
    #if len(variants) > 0:
    #    transcriptInsertions = createIndelBedtool(transcripts, "I", maxLen)
    #    varDict = intersectWithVariants(transcriptInsertions, variants)

    for t in transcripts.keys():
        t = transcripts[t]

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
    
        # Iterate over operations to sequence and repair mismatches and microindels
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
                    # Variant-aware mode
                    if len(variants) > 0:

                        # Check if the insertion is in the optional variant catalog.
                        if ID in variants:
                            # The insertion perfectly matches a variant. Do not correct it.
                            comment = "DidNotCorrect_insertion_at_" + currPos + "_becauseVariantMatch"
                            Transcript2.updateLog(t, comment)
                            continue

                    else: # Correct insertion
                        # Add comment to log indicating that we made a change to the transcript
                        comment = "Corrected_deletion_at_" + ID
                        Transcript2.updateLog(t, comment)
                        # Subtract the inserted bases by skipping them. GenomePos stays the same, as does MVal
                        seqPos += ct
                else: # Move on without correcting insertion because it is too big
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

            # N, H, and D operations are cases where the transcript sequence is missing bases that are in the reference genome.
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
        #t.NM, t.MD = t.getNMandMDFlags(genome) # May not be necessary since MD tag lacks I
    return    

def correctDeletions(transcripts, genome, variants, maxLen):
    # This function corrects deletions up to size maxLen using the reference genome. If a variant file was provided, correction will be SNP-aware.

    # If variants are provided, get all transcript insertions and intersect them with the variant catalog
    #if len(variants) > 0:
    #    transcriptDeletions = createIndelBedtool(transcripts, "D", maxLen)
    #    varDict = intersectWithVariants(transcriptDeletions, variants)

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
                if ct <= maxLen:
                    ID = "_".join([t.CHROM, str(genomePos), str(genomePos + ct - 1)])
                    # Variant-aware mode
                    if len(variants) > 0:

                        # Check if the deletion is in the optional variant catalog.
                        if ID in variants:
                            # The deletion perfectly matches a variant. Do not correct it.
                            comment = "DidNotCorrect_deletion_at_" + currPos + "_becauseVariantMatch"
                            Transcript2.updateLog(t, comment)
                            continue
                    else: # Correct deletion if we're not in variant-aware mode
                        # Add comment to log indicating that we made a change to the transcript
                        comment = "Corrected_deletion_at_" + ID
                        Transcript2.updateLog(t, comment)
                        # Add the missing reference bases
                        refBases = genome.sequence({'chr': t.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                        newSeq = newSeq + refBases
                        genomePos += ct
                        MVal += ct
                # Deletion is too big to fix
                else:
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

def createIndelBedtool(transcripts, operation, maxLen):
    # This function iterates through the transcripts and creates a BedTool object containing each insertion or deletion (depending on whether "I" or "D" is specified in the operation variable).
    # The purpose is to create a data structure that can be intersected with a VCF file of variants.

    o = open("tmp_indel.bed", 'w')
    for tID in transcripts.keys():
        t = transcripts[tID]
        cigarOps,cigarCounts = t.splitCIGAR()
    
        # Check for operation. If none are present, we can skip this transcript
        if operation not in t.CIGAR : continue

        # Start at position in the genome where the transcript starts.
        genomePos = t.POS - 1 # Make it zero-based since we're making a bed file

        # Iterate over CIGAR operations to find positions of the operations that we care about
        for op,ct in zip(cigarOps, cigarCounts):

            if op == "D":
                if ct <= maxLen and operation == "D":
                    bedChrom = t.CHROM
                    bedStart = genomePos
                    bedEnd = bedStart + ct
                    bedName = tID
                    bedStrand = t.strand
                    o.write("\t".join([bedChrom, str(bedStart), str(bedEnd), bedName, ".", bedStrand]) + "\n")
                genomePos += ct
            if op == "I":
                if ct <= maxLen and operation == "I":
                    bedChrom = t.CHROM
                    bedStart = genomePos
                    bedEnd = bedStart + ct
                    bedName = tID
                    bedStrand = t.strand
                    o.write("\t".join([bedChrom, str(bedStart), str(bedEnd), bedName, ".", bedStrand]) + "\n")
            if op in ["M", "N", "H"]:
                genomePos += ct

    o.close()
    os.system('sort -k1,1 -k2,2n tmp_indel.bed > sorted_tmp_indel.bed')
    indel = pybedtools.BedTool("sorted_tmp_indel.bed")
    os.system('rm tmp_indel.bed')
    #os.system('rm sorted_tmp_indel.bed')

    return indel

def correctMismatches(transcripts, genome, variants):
    # This function corrects mismatches in the sequences. If a variant file was provided, correction will be SNP-aware.
   
    for t in transcripts.keys():
        t = transcripts[t]

        origSeq = t.SEQ
        origCIGAR = t.CIGAR
        origMD = t.MD

        # Check for mismatches. If none are present, we can skip this transcript
        if any(i in origMD.upper() for i in 'ACTGN') == False : continue

        newCIGAR = ""
        newSeq = ""
        MVal = 0
        seqPos = 0
        genomePos = t.POS

        # Merge CIGAR and MD tag information so that we have the locations of all insertions, deletions, and mismatches
        mergeOperations, mergeCounts = t.mergeMDwithCIGAR()

        # Iterate over operations to sequence and repair mismatches and microindels
        for op,ct in zip(mergeOperations, mergeCounts):
 
            if op == "M":
                 newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                 MVal += ct
                 seqPos += ct
                 genomePos += ct

            # This denotes a mismatch
            if op == "X":
                # Check if the position matches a variant in the optional variant catalog. If yes, don't try to fix it.
                isSNP = False
                try:
                    currBase = origSeq[seqPos]
                    if currBase in variants[(t.CHROM + "_" + str(genomePos))]: 
                        # Add comment to log indicating that we declined to change the transcript
                        comment = "DidNotCorrect_mismatch_at_" + t.CHROM + ":" + str(genomePos)
                        Transcript2.updateLog(t, comment)
                        # Keep the base as-is
                        isSNP = True
                        newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                        MVal += ct
                        seqPos += ct
                        genomePos += ct
                    else: isSNP = False
                except: pass
                if isSNP == False:
                    # Add comment to log indicating that we made a change to the transcript
                    comment = "Corrected_mismatch_at_" + t.CHROM + ":" + str(genomePos)
                    Transcript2.updateLog(t, comment)

                    # Change sequence base to the reference base at this position
                    newSeq = newSeq + genome.sequence({'chr': t.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
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

        t.CIGAR = newCIGAR
        t.SEQ = newSeq
        t.NM, t.MD = t.getNMandMDFlags(genome)
    return 

def correctMismatchesAndIndels(transcripts, genome, variants, maxLen):
    # This function corrects mismatches and indels up to size maxLen using the reference genome. If a variant file was provided, correction will be SNP-aware.
    
    for t in transcripts.keys():
        t = transcripts[t]

        origSeq = t.SEQ
        origCIGAR = t.CIGAR
        origMD = t.MD
        genomePos = t.POS

        print origCIGAR
        print origMD

        # Check for deletions, insertions, and mismatches. If none are present, we can skip this transcript
        if ("D" not in origCIGAR) and ("I" not in origCIGAR) and any(i in origMD.upper() for i in 'ACTGN') == False : continue
      
        newCIGAR = ""
        newSeq = ""
        MVal = 0
        seqPos = 0
        genomePos = t.POS
       
        # Merge CIGAR and MD tag information so that we have the locations of all insertions, deletions, and mismatches
        mergeOperations, mergeCounts = t.mergeMDwithCIGAR()
        
        # Iterate over operations to sequence and repair mismatches and microindels
        for op,ct in zip(mergeOperations, mergeCounts):
            if op == "M":
                 newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                 MVal += ct
                 seqPos += ct
                 genomePos += ct
            # This denotes a mismatch
            if op == "X": 
                # TODO: check if the deletion is in the optional variant catalog. If yes, don't try to fix it.
                # Change sequence base to the reference base at this position
                newSeq = newSeq + genome.sequence({'chr': t.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                seqPos += ct # skip the original sequence base
                genomePos += ct # advance the genome position
                MVal += ct
               
            if op == "D":
                 # TODO: check if the deletion is in the optional variant catalog. If yes, don't try to fix it.
                if ct <= maxLen:
                    # Add the missing reference bases
                    refBases = genome.sequence({'chr': t.CHROM, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                    newSeq = newSeq + refBases
                    genomePos += ct
                    MVal += ct
                else:
                    # End any ongoing match
                    MVal, newCIGAR = endMatch(MVal, newCIGAR)
                    newCIGAR = newCIGAR + str(ct) + op
                    genomePos += ct
            if op == "I":
            # TODO: check if the insertion is in the optional variant catalog. If yes, don't try to fix it.
                if ct <= maxLen:
                    # Subtract the inserted bases by skipping them. GenomePos stays the same, as does MVal
                    seqPos += ct
                else:
                    # End any ongoing match
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
                        
def endMatch(MVal, newCIGAR):
    # Function to end an ongoing match during error correction, updating new CIGAR and MD tag strings as needed.
    # MVal keeps track of how many matched bases we have at the moment so that when errors are fixed, adjoining matches can be combined.
    # When an intron or unfixable error is encountered, it is necessary to end the match by outputting the current MVal and then setting the variable to zero.

    if MVal > 0:
        newCIGAR = newCIGAR + str(MVal) + "M"
        MVal = 0

    return MVal, newCIGAR

def cleanNoncanonical(transcripts, annotatedJunctions, genome, n):
    # Iterate over noncanonical transcripts. Determine whether each end is within n basepairs of an annotated junction.
    # If it is, run the rescue function on it. If not, discard the transcript.

    o = open("tmp_nc.bed", 'w')
    salvageableNCJns = 0
    totNC = len(transcripts)
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
    os.system('sort -k1,1 -k2,2n tmp_nc.bed > sorted_tmp_nc.bed')
    nc = pybedtools.BedTool("sorted_tmp_nc.bed")
    jnMatches = str(nc.closest(annotatedJunctions, s=True, D="ref", t="first")).split("\n")

    os.system("rm tmp_nc.bed")
    os.system("rm sorted_tmp_nc.bed")

    # Iterate over splice junction boundaries and their closest canonical match. 
    for match in jnMatches:
        if len(match) == 0: continue
        match = match.split('\t')
        d = int(match[-1])
        transcriptID, spliceJnNum, side = match[3].split("__")

        # Only attempt to rescue junction boundaries that are within n bp of an annotated junction
        if abs(d) > n:
            continue

        currTranscript = transcripts[transcriptID]
        currJunction = currTranscript.spliceJunctions[int(spliceJnNum)]
        currIntronBound = currJunction.bounds[int(side)]
        rescueNoncanonicalJunction(currTranscript, currJunction, currIntronBound, d, genome)

        # Add comment to log indicating that we made a change to the transcript
        comment = "Corrected_SpliceJn_" + str(spliceJnNum) + ":" + str(side) + "_" + str(n) + "bp"
        Transcript2.updateLog(currTranscript, comment)

        #t.NM, t.MD = t.getNMandMDFlags(genome)
    return

def rescueNoncanonicalJunction(transcript, spliceJn, intronBound, d, genome):

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
            spliceJn.end = intronBound.pos

        if d < 0: # Need to subtract from end of exon sequence. Case 3
            exonSeqs[targetExon] = exon[0:d]
            intronBound.pos += d
            spliceJn.start = intronBound.pos
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
            spliceJn.end = intronBound.pos
        if d > 0: # Need to subtract from start of exon sequence. Case 4
            exonSeqs[targetExon] = exon[d:]
            intronBound.pos += d
            spliceJn.end = intronBound.pos
        # Modify exon string
        exonCIGARs[targetExon] = editExonCIGAR(exonCIGARs[targetExon], 0, -d)
        intronCIGARs[targetJn] += d
    transcript.SEQ = ''.join(exonSeqs)

    # Paste together the new CIGAR string
    newCIGAR = ""
    for i in range(0,len(intronCIGARs)):
        newCIGAR = newCIGAR + exonCIGARs[i] + str(intronCIGARs[i]) + "N"
    newCIGAR = newCIGAR + exonCIGARs[-1]
    transcript.CIGAR = newCIGAR    
    intronBound.isCanonical = True
    return

 
def splitCIGAR(CIGAR):
    # Takes CIGAR string from SAM and splits it into two lists: one with capital letters (match operators), and one with the number of bases

    matchTypes = re.sub('[0-9]', " ", CIGAR).split()
    matchCounts = re.sub('[A-Z]', " ", CIGAR).split()
    matchCounts = [int(i) for i in matchCounts]

    return matchTypes, matchCounts

def editExonCIGAR(exon, side, nBases):
    # Given an exon CIGAR string, this function adds or subtracts bases from the start (side = 0) or end (side = -1)
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
        



main()
