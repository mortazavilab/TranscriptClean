# TranscriptClean
TranscriptClean is a Python program that corrects mismatches, microindels, and noncanonical splice junctions in long reads that have been mapped to the genome. It is designed for use with sam files from the PacBio Iso-seq and Oxford Nanopore transcriptome sequencing technologies. A variant-aware mode is available for users who want to avoid correcting away known variants in their data.

Note: At the present time, TranscriptClean does not work on SAM files that use X/= operators rather than M to represent matches in the CIGAR field. We are working on adding support for this in a future version.

## Installation

To install, clone this repo and install using pip:
```bash
git clone git@github.com:mortazavilab/TranscriptClean.git
cd TranscriptClean
pip install -e .
```

<!-- The current TranscriptClean version is designed to be run with Python >= 3.7. It requires Bedtools to be installed, as well as Python modules pybedtools and pyfasta. These can be found at the links listed below:

```
conda install -c bioconda python=3.7 pyranges samtools pyfaidx
```

* pyfaidx (v0.7.1): https://pypi.org/project/pyfaidx/
* Samtools (v1.9): https://github.com/samtools/samtools/releases/
* PyRanges: https://pyranges.readthedocs.io
* pybedtools (optional for variants) (v0.7.8): https://daler.github.io/pybedtools/
* Bedtools (optional for variants) (v2.25.0): http://bedtools.readthedocs.io/en/latest/content/installation.html -->

In addition, R (tested with v.3.3.2) is needed to run the visualization script, generate_report.R.

<!-- To install TranscriptClean, simply download the files using Github's "Download ZIP" button, then unzip them in the directory where you would like to install the program. Alternately, you can download a specific version of the program from the Releases tab. The TranscriptClean script can now be run directly from the command line- just include the path. -->

## Usage
TranscriptClean is run from the command line as follows. Please note that releases 2.0+ can be run in multithreaded fashion. For fastest performance, we recommend sorting your input SAM file.
**For additional details and examples, please see the Wiki section.**

`transcriptclean --sam transcripts.sam --genome hg38.fa --outprefix /my/path/outfile`


### Basic Options
| Option            | Shortcut  | Description
|------------------ | --------- | -----------------------------------------------------------------------------------------------------------------------
| --help            | -h        | Print a list of the input options with descriptions
| --sam 	    | -s        | Input sam file (mandatory). The aligner used to create it must be splice aware if you want to correct splice junctions.
| --genome          | -g        | Reference genome fasta file (mandatory). Should be the same one used during alignment to generate the sam file.
| --threads         | -t        | Number of threads to run program with. Default = 1.
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
| --primaryOnly       | n/a       | If this option is set, TranscriptClean will only output primary mappings of transcripts (ie it will filter out unmapped and multimapped lines from the SAM input.)
| --canonOnly       | n/a       | If this option is set, TranscriptClean will only output transcripts that are either canonical or that contain annotated noncanonical junctions to the clean SAM and Fasta files at the end of the run.

### Other options that may help tune performance
| Option              | Shortcut  | Description
|-------------------- | --------- | ---------------------------------------------------------------------------------------------------------------------
| --tmpDir           | n/a       | If you would like the tmp files to be written somewhere different than the final output, provide the path to that location here. For example, a tmp directory on the local drive of a compute node.
| --bufferSize        | n/a       | Number of lines to output to file at once by each thread during run. Default = 100
| --deleteTmp         | n/a       | If this option is set, the temporary directory generated by TranscriptClean (TC_tmp) will be removed at the end of the run.


## Output files
TranscriptClean outputs the following files:
* SAM file of corrected transcripts. Unmapped/non-primary transcript alignments from the input file are included in their original form.
* Fasta file of corrected transcript sequences. Unmapped transcripts from the input file are included in their original form.
* Transcript error log file (.TE.log): Each row represents a potential error in a given transcript. The column values track whether the error was corrected or not and why.
* Transcript log file (.log): Each row represents a transcript. The columns track the mapping status of the transcript, as well as how many errors of each type were found and corrected/not corrected in the transcript.

## Credit
Please cite our paper when using TranscriptClean:

Dana Wyman, Ali Mortazavi, TranscriptClean: variant-aware correction of indels, mismatches and splice junctions in long-read transcripts, Bioinformatics, Volume 35, Issue 2, 15 January 2019, Pages 340–342, https://doi.org/10.1093/bioinformatics/bty483
