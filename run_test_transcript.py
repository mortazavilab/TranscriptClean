# THis testing script takes in a sam file and tries processing it as a way to make sure that the transcript.py class is working as expected.

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
    parser.add_option("--n", dest = "nLines", help = "Number of entries to print out from the test file. -1 will print them all.",
                      type = "int", default = "-1")
    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    nLines = int(options.nLines)
    count = 0

    # Read in the reference genome. Treat coordinates as 0-based 
    print "Reading genome .............................."
    genome = Fasta(options.refGenome)
   
    print "Processing SAM file ........................."
    with open(options.sam, 'r') as f:
        for line in f:
            if nLines > -1 and count >= nLines:
                break

            line = line.strip()
            if line.startswith("@"):
                continue

            print "---------------------------------" + "\n" + line + "\n"
            t = Transcript(line, genome)
            if int(t.FLAG) > 16:
            #    count +=1
                continue
            # Skip unmapped transcripts altogether
            if t.CHROM == "*":
            #    count += 1
                continue

            print Transcript.printableSAM(t, genome) + "\n"
            count += 1

main()
