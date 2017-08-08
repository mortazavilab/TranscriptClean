
from transcript import Transcript
from spliceJunction import SpliceJunction
from optparse import OptionParser

def getOptions():
    parser = OptionParser()
    parser.add_option("--f", dest = "infile", help = "Input file",
                      metavar = "FILE", type = "string", default = "")
    parser.add_option("--o", dest = "outfile",
                      help = "output file", metavar = "FILE", type = "string", default = "out")
    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    
    header = ""
    with open(options.infile, 'r') as f:
        for line in f:
	    line = line.strip()
            if line.startswith("@"):
                header = header + line + "\n"
                continue
            t = Transcript(line)
            t.printSpliceJunctions()
            exit()
main()
