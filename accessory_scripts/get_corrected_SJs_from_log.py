# This script parses the TE log from TranscriptClean to extract the formely noncanonical 
# splice junctions in BED format.

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

o = open(outfile, 'w')

with open(infile, 'r') as f:
    for line in f:
        if line.startswith("Transcript"): continue

        info = line.strip().split("\t")
        if info[2] == "NC_SJ_boundary" and info[4] == "Corrected": 
            print info

            dist_0, dist_1 = info[3].split("_")
            dist_0 = int(dist_0)
            dist_1 = int(dist_1)
  
            positionInfo = info[1].split("_")
            chrom = "_".join(positionInfo[0::-2])
            start = int(positionInfo[-2]) + dist_0 - 1
            end = int(positionInfo[-1]) + dist_1

            outstr = "\t".join([chrom, str(start), str(end)]) + "\n"
            o.write(outstr)

o.close()
