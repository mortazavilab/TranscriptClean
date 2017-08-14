
from transcript import Transcript
from spliceJunction import SpliceJunction
from optparse import OptionParser
import pybedtools
from pyfasta import Fasta
import os

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
    cleanNoncanonical(noncanTranscripts, annotatedSpliceJns)

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

def cleanNoncanonical(transcripts, annotatedJunctions):
    # Iterate over noncanonical transcripts. Determine whether each end is within 5 basepairs of an annotated junction.
    # If it is, run the rescue function on it. If not, discard the transcript.

    o = open("tmp_nc.bed", 'w')
    #noncanonicalJns = 0
    salvageableNCJns = 0
    totNC = len(transcripts)
    for tID in transcripts.keys():
        t = transcripts[tID]
        jns = t.spliceJunctions
        for jn in jns:
            if jn.isCanonical == True:
                continue
            
            #noncanonicalJns += 1

            # Get BedTool object for start of junction
            start = SpliceJunction.getBED(jn, "start")
            #d1 = closestAnnotatedJn(start, annotatedJunctions)
            o.write(start + "\n")            

            # Get BedTool object for start of junction
            end = SpliceJunction.getBED(jn, "end")
            #d2 = closestAnnotatedJn(end, annotatedJunctions)
            o.write(end + "\n")

       #     if abs(d1) <= 5 and abs(d2) <=5:
       #         salvageableNCJns += 1
       #         SpliceJunction.changeToCanonical(jn)
       # if Transcript.recheckCanonical(t) == True:
       #     cleanTranscripts.append(t)

    o.close()
    os.system('sort -k1,1 -k2,2n tmp_nc.bed > sorted_tmp_nc.bed')
    nc = pybedtools.BedTool("sorted_tmp_nc.bed")
    matches = str(nc.closest(annotatedJunctions, s=True, D="ref", t="first")).split("\n")
    
    for entry in matches:
        if len(entry) == 0: continue
        entry = entry.split('\t')
        d = int(entry[-1])
        transcriptID, spliceJnNum, side = entry[3].split("__")
        if abs(d) > 5:
            transcripts.pop(transcriptID, None)
        

    #print "Total noncanonical junctions: " + str(noncanonicalJns)
    #percentSalvageableJn = round(float(salvageableNCJns)*100/noncanonicalJns, 2)
    #print "Number of salvageable noncanonical junctions: " + str(salvageableNCJns) + " (" + str(percentSalvageableJn) + "%)"
    percentSalvageableT = round(float(len(transcripts))*100/totNC, 2)
    print "Number of salvageable noncanonical transcripts: " + str(len(transcripts)) + " (" + str(percentSalvageableT) + "%)"

def closestAnnotatedJn(jn, annotatedJunctions):
    # This function accepts a noncanonical splice junction BedTool object and finds the closest annotated junction. It also
    # reports the distance from said junction.
    
    closest = str(jn.closest(annotatedJunctions, s=True, D="ref", t="first" )[0]).strip().split()
    return int(closest[-1])
    
    

def rescueNoncanonicalJunction(jn, d):
    # This function converts a noncanonical splice junction to a canonical junction that is <= 5 bp away.
    # To do this, it is necessary to
    # (1) Change the sam sequence to the reference
    # (2) Potentially change the mapping quality? (Not sure how yet)
    # (3) Change the CIGAR string
    # (4) Change the splice junction introns
    # (5) Change the splice junction string   
    pass
main()
