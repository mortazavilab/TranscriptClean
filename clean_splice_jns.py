
from transcript import Transcript
from spliceJunction import SpliceJunction
from intronBound import IntronBound
from optparse import OptionParser
import pybedtools
from pyfasta import Fasta
import os
import re

def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "sam", help = "Input file",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--g", dest = "refGenome", help = "Reference genome fasta file. Should be the same one used to generate the sam file.",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--s", dest = "spliceAnnot", help = "Splice junction file obtained by mapping Illumina reads to the genome using STAR. More formats may be supported in the future.",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--o", dest = "outfile",
                      help = "output file", metavar = "FILE", type = "string", default = "out")
    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
   
    # Read in the reference genome. Treat coordinates as 0-based 
    print "Reading genome .............................."
    genome = Fasta(options.refGenome)
    print "Processing SAM file ........................."
    header, canTranscripts, noncanTranscripts = processSAM(options.sam, genome)
    print "Processing annotated splice junctions ........"
    annotatedSpliceJns = processSpliceAnnotation(options.spliceAnnot)
    

    totalT = len(canTranscripts) +  len(noncanTranscripts)
    print "Total transcripts: " + str(totalT)
    print "Total noncanonical transcripts: " + str(len(noncanTranscripts))
    cleanNoncanonical(noncanTranscripts, annotatedSpliceJns, genome)

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
            t = Transcript(line, genome)
            if t.isCanonical == True:
                canTranscripts[t.QNAME] = t
            else:
                noncanTranscripts[t.QNAME] = t
    return header, canTranscripts, noncanTranscripts

def processSpliceAnnotation(annotFile):
    # This function reads in the tab-separated STAR splice junction file and creates a bedtools object
    
    bedstr = ""
    o = open("tmp.bed", 'w')
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
    os.system('sort -k1,1 -k2,2n tmp.bed > tmp2.bed')
    bt = pybedtools.BedTool("tmp2.bed")
    return bt    

def cleanNoncanonical(transcripts, annotatedJunctions, genome):
    # Iterate over noncanonical transcripts. Determine whether each end is within 5 basepairs of an annotated junction.
    # If it is, run the rescue function on it. If not, discard the transcript.

    o = open("tmp_nc.bed", 'w')
    salvageableNCJns = 0
    totNC = len(transcripts)
    for tID in transcripts.keys():
        t = transcripts[tID]
        bounds = Transcript.getAllIntronBounds(t)
        
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
    os.system("rm tmp.bed")
    os.system("rm tmp2.bed")

    # Iterate over splice junction boundaries and their closest canonical match. 
    for match in jnMatches:
        if len(match) == 0: continue
        match = match.split('\t')
        d = int(match[-1])
        transcriptID, spliceJnNum, side = match[3].split("__")
        
        # Only attempt to rescue junction boundaries that are within 5 bp of an annotated junction
        if abs(d) > 5:
            transcripts.pop(transcriptID, None)
            continue
        
        currTranscript = transcripts[transcriptID]
        currJunction = currTranscript.spliceJunctions[int(spliceJnNum)]
        currIntronBound = currJunction.bounds[int(side)]
        rescueNoncanonicalJunction(currTranscript, currJunction, currIntronBound, d)

    print transcripts["PB.9671.7|chr20:46118311-46129508(+)|c38618/f3p1/1373"].CIGAR

def rescueNoncanonicalJunction(transcript, spliceJn, intronBound, d):
    # This function converts a noncanonical splice junction to a canonical junction that is <= 5 bp away.
    # To do this, it is necessary to
    # (1) Edit the sam sequence using the reference
    # (2) Potentially change the mapping quality? (Not sure how yet)
    # (3) Change the CIGAR string
    # (4) Change the splice junction intron coordinates
    # (5) Change the splice junction string   
    
    seq = transcript.SEQ
  
    correctCIGAR(transcript, spliceJn.jnNumber, intronBound.bound, d)
    #print transcript.QNAME
    #print spliceJn.jnNumber
    #print intronBound
    #print d

def correctCIGAR(transcript, targetIntron, intronBound, d):
    # This function modifies the CIGAR string of a transcript to make it reflect the conversion of a particular 
    # splice junction from noncanonical to canonical

    CIGAR = transcript.CIGAR

    # The exon we need to modify depends on which end of the splice junction we are rescuing
    # If we are rescuing the left side, we will modify the exon with the same number as the intron
    # If we are rescuing the right, we will modify exon number intron+1
    targetExon = targetIntron + intronBound

    
    matchTypes, matchCounts = splitCIGAR(CIGAR)
    newCIGAR = ""
    currIntron = 0
    currExon = 0
    for operation,c in zip(matchTypes, matchCounts):
        if operation == "M":   
            if currExon == targetExon:
                if (intronBound == 0 and d > 0) or (intronBound == 1 and d < 0):
                    # Under these conditions, the exon length will increase. Adjust the match count accordingly
                    c = c + abs(d)
                if (intronBound == 0 and d < 0) or (intronBound == 1 and d > 0):
                    # Under these conditions, the exon length will decrease. 
                    c = c - abs(d)
        if operation == "N": 
            if currIntron == targetIntron:
                if (intronBound == 0 and d > 0) or (intronBound == 1 and d < 0):
                    # Under these conditions, the intron length will decrease. Adjust the match count accordingly
                    c = c - abs(d)
                if (intronBound == 0 and d < 0) or (intronBound == 1 and d > 0):
                    # Under these conditions, the intron length will increase. 
                    c = c + abs(d)
            currIntron += 1
            currExon += 1 
        newCIGAR = newCIGAR + str(c) + operation
    # Update the transcript
    transcript.CIGAR = newCIGAR
    return 

def splitCIGAR(CIGAR):
    # Takes CIGAR string from SAM and splits it into two lists: one with capital letters (match operators), and one with the number of bases

    matchTypes = re.sub('[0-9]', " ", CIGAR).split()
    matchCounts = re.sub('[A-Z]', " ", CIGAR).split()
    matchCounts = [int(i) for i in matchCounts]

    return matchTypes, matchCounts

def getExonSeqs(seq, CIGAR):
    # Given a sequence (from a sam entry) and a CIGAR string, this function splits the sequence into substrings 
    # corresponding to the size and order of the CIGAR string M operations. It returns only the exon (M) sequences
    # Example: seq = ATCGATTCA and CIGAR 2M4N3M:
    #          Yields matchTypes = [ M, N, M ] and matchCounts = [ 2, 4, 3 ] 
    #          Result = [ AT, TCA ]

    matchTypes, matchCounts = splitCIGAR(CIGAR)
    result = []
   
    curr = 0
    for m,c in zip(matchTypes, matchCounts):
        if m != "M": continue
        subSeq = seq[curr:curr + c]
        curr = curr + c
        result.append(subSeq)
    return result


main()
