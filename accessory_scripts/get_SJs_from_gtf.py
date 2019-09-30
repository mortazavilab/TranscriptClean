# This script reads in a GTF transcript annotation and extracts the splice 
# junctions. Exons must be in order.
# The output format is designed to match the STAR SJ file output format

from optparse import OptionParser
from pyfasta import Fasta

def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "infile", help = "Input GTF file",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--g", dest = "genomeFile", help = "Reference genome",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--minIntronSize", "-m", dest = "minIntron", default = 21,
                      help = "Minimum size of intron to consider a junction. Default: 21 bp.")
    parser.add_option("--o", dest = "outfile",
                      help = "output file", metavar = "FILE", type = "string", default = "out.txt")
    (options, args) = parser.parse_args()
    return options

def formatSJOutput(currExon, prev_exonEnd, genome, minIntron):
    chromosome = currExon[0]
    strand = currExon[6]
    if strand == "+": 
        strand = "1"
        intron_start = int(prev_exonEnd) + 1
        intron_end = int(currExon[3]) - 1
        intronMotif = getIntronMotif(chromosome, intron_start, intron_end, genome)
    elif strand == "-":
        strand = "2"
        intron_start = int(currExon[4]) + 1  #int(currExon[3]) + 1
        intron_end = int(prev_exonEnd) - 1    #int(prev_exonEnd) - 1
        intronMotif = getIntronMotif(chromosome, intron_start, intron_end, genome)
    if abs(intron_end - intron_start + 1) < minIntron:
        return None 
    intronMotif = getIntronMotif(chromosome, intron_start, intron_end, genome)
    annotationStatus = "1"
    nUniqueReads = "NA"
    nMultiReads = "NA"
    maxSpliceOverhang = "NA"
    
    return "\t".join([chromosome, str(intron_start), str(intron_end), strand, intronMotif, annotationStatus, nUniqueReads, nMultiReads, maxSpliceOverhang])

def getIntronMotif(chrom, start, end, genome):
   startBases = genome.sequence({'chr': chrom, 'start': start, 'stop': start + 1}, one_based=True)
   endBases = genome.sequence({'chr': chrom, 'start': end - 1, 'stop': end}, one_based=True)
   motif = (startBases + endBases).upper() 

   if motif == "GTAG":
       return "21"
   elif motif == "CTAC":
       return "22"
   elif motif == "GCAG":
       return "23"
   elif motif == "CTGC":
       return "24"
   elif motif == "ATAC":
       return "25"
   elif motif == "GTAT":
       return "26"
   else:
       return "20" 

if __name__ == "__main__":

    junctions_seen = {}

    # Read input arguments
    options = getOptions()
    gtf = options.infile
    genome = Fasta(options.genomeFile) 
    minIntron = int(options.minIntron)
    o = open(options.outfile, 'w')

    # Read in the GTF
    prev_transcriptID = ""
    prev_exonEnd = 0
    with open(gtf, 'r') as f:
        for line in f:

            # Prep
            line = line.strip()
            
            # Ignore header
            if line.startswith("#"):
                continue 
         
            # Split GTF line on tab
            info = line.split("\t")

            # Ignore entries that are not exons
            if info[2] != "exon":
                continue

            # Extract transcriptID and exonID from description field
            description = info[-1]

            # Skip entries that lack a transcript ID
            if "transcript_id" not in description:
                continue
            
            transcriptID = (description.split("transcript_id ")[1]).split('"')[1]
            strand = info[6]

            if transcriptID != prev_transcriptID:
                # Start new transcript
                prev_transcriptID = transcriptID
                if strand == "+":
                    prev_exonEnd = info[4]
                else:
                    prev_exonEnd = info[3]
            else: 
                # Output the current junction 
                spliceJn = formatSJOutput(info, prev_exonEnd, genome, minIntron)
                if strand == "+":
                    prev_exonEnd = info[4]
                else:
                    prev_exonEnd = info[3]
                if spliceJn != None:
                    if spliceJn not in junctions_seen:
                        o.write(spliceJn + "\n")
                        junctions_seen[spliceJn] = 1
    o.close()
            
