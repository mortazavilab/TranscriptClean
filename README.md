# TranscriptClean
TranscriptClean is a Python program that corrects mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome. It is designed for use with sam files from the PacBio Iso-seq and Oxford Nanopore transcriptome sequencing technologies. A variant-aware mode is available for users who want to avoid correcting away known variants in their data.

Note: At the present time, TranscriptClean does not work on SAM files that use X operators rather than M to represent matches in the CIGAR field. We are working on adding support for this in a future version. 

## Installation
TranscriptClean is designed to be run with Python version 2.7. It requires Bedtools to be installed, as well as Python modules pybedtools and pyfasta. These can be found at the links listed below:
* Bedtools (v2.25.0): http://bedtools.readthedocs.io/en/latest/content/installation.html
* pybedtools (v0.7.8): https://daler.github.io/pybedtools/
* pyfasta (v0.5.2): https://pypi.python.org/pypi/pyfasta/

In addition, R (v.3.3.2 recommended) is needed to run the visualization script, generate_report.R.

To install TranscriptClean, simply download the files using Github's "Download ZIP" button, then unzip them in the directory where you would like to install the program. Alternately, you can download a specific version of the program from the Releases tab. The TranscriptClean script can now be run directly from the command line- just include the path. 

## Usage 
TranscriptClean is run from the command line as follows. For additional details and examples, please see the Wiki section.

`python TranscriptClean.py --sam transcripts.sam --genome hg38.fa --outprefix /my/path/outfile`


### Basic Options
| Option            | Shortcut  | Description
|------------------ | --------- | ----------------------------------------------------------------------------------------------------------------------- 
| --help            | -h        | Print a list of the input options with descriptions
| --sam 	    | -s        | Input sam file (mandatory). The aligner used to create it must be splice aware if you want to correct splice junctions.
| --genome          | -g        | Reference genome fasta file (mandatory). Should be the same one used during alignment to generate the sam file.
| --outprefix       | -o        | Prefix for the output files. Default = "out".

### Options that control run mode
| Option              | Shortcut  | Description
|-------------------- | --------- | ---------------------------------------------------------------------------------------------------------------------
| --dryRun            | n/a       | Include this option to run an inventory of all indels in the data without performing any correction. Useful for selecting maxLenIndel and maxSJOffset size
| --correctMismatches | -m        | If set to false, TranscriptClean will skip mismatch correction. Default = True.
| --correctIndels     | -i        | If set to false, TranscriptClean will skip indel correction. Default = True.
| --variants          | -v        | Optional: VCF-formatted file of variants to avoid correcting (this enables variant-aware correction). Irrelevant if correctMismatches is set to false. 
| --spliceJns         | -j        | High-confidence splice junction file obtained by mapping Illumina short reads to the genome using STAR. More formats may be supported in the future. This file is necessary if you want to correct noncanonical splice junctions.
| --maxLenIndel       | n/a       | Maximum size indel to correct. Default = 5 bp.
| --maxSJOffset       | n/a       | Maximum distance from annotated splice junction to correct. Default = 5 bp.         

## Output files
TranscriptClean outputs the following files:
* SAM file of corrected transcripts. Unmapped/non-primary transcript alignments from the input file are included in their original form.
* Fasta file of corrected transcript sequences. Unmapped transcripts from the input file are included in their original form.
* Transcript error log file (.TE.log): Each row represents a potential error in a given transcript. The column values track whether the error was corrected or not and why.
* Transcript log file (.log): Each row represents a transcript. The columns track the mapping status of the transcript, as well as how many errors of each type were found and corrected/not corrected in the transcript.

## Credit
Please cite our paper when using TranscriptClean:

Dana Wyman, Ali Mortazavi. TranscriptClean: variant-aware correction of indels, mismatches and splice junctions in long-read transcripts, Bioinformatics, 15 June 2018, https://doi.org/10.1093/bioinformatics/bty483


