# TranscriptClean
TranscriptClean is a Python program that corrects mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome. It is designed for use with sam files from the PacBio Iso-seq and Oxford Nanopore transcriptome sequencing technologies. A variant-aware mode is available for users who want to avoid correcting away SNPs in their data.

## Installation
TranscriptClean is designed to be run with Python version 2.7. In addition, TranscriptClean requires Bedtools to be installed, as well as Python modules pybedtools, pyfasta, and re (regular expressions). These can be found at the links listed below:
* Bedtools: http://bedtools.readthedocs.io/en/latest/content/installation.html
* pybedtools: https://daler.github.io/pybedtools/
* pyfasta: 
* re: https://pypi.python.org/pypi/re2/


## Usage 
TranscriptClean is run from the command line as follows:

`python TranscriptClean.py --sam transcripts.sam --genome hg38.fa --outprefix /my/path/outfile`

Option | Shortcut | Description
-------------- | ------------- | ------------- 
--help 	| -h | Input sam file (mandatory). The aligner used to create it must be splice aware if you want to correct splice junctions.
--genome | -g | Reference genome fasta file (mandatory). Should be the same one used during alignment to generate the sam file.
--outprefix | -o | Prefix for the output files. Default = "out".


