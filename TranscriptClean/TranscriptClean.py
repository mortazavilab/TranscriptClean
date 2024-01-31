# TranscriptClean
# Author: Dana Wyman
# 3/1/2018
# -----------------------------------------------------------------------------
# Mismatches and microindels in long reads are corrected in a variant-aware
# fashion using the reference genome and a VCF file of whitelisted variants.
# Noncanonical splice junctions can also be corrected using a file of reference
# splice sites.

# TC Classes
from TranscriptClean.transcript import Transcript
from TranscriptClean.transcript import check_seq_and_cigar_length
from TranscriptClean.spliceJunction import *
from TranscriptClean.intronBound import IntronBound
from optparse import OptionParser
from TranscriptClean.dstruct import Struct

# Other modules
from pyfaidx import Fasta
import os
import multiprocessing as mp
from math import ceil
import re
from copy import copy
import warnings

# Runtime profiling
import cProfile
import time
from datetime import timedelta


def main():
    options = getOptions()
    sam_file = options.sam
    n_threads = options.n_threads
    header, sam_chroms, sam_chunks = split_SAM(sam_file, n_threads)
    validate_chroms(options.refGenome, options.variantFile, sam_chroms)

    # Run the processes. Outfiles are created within each process
    processes = []

    # run_chunk(sam_chunks[0], options, header)

    # for i in range(n_threads):
    # only run in range n_threads in case multiproc allocation
    # doesn't need all threads
    for i in range(len(sam_chunks)):
        if options.dryRun == True:
            t = mp.Process(target=run_chunk_dryRun,
                           args=(sam_chunks[i], options))
        else:
            t = mp.Process(target=run_chunk, args=(
                sam_chunks[i], options, header))

        # Run the process
        processes.append(t)
        t.start()

    # Wait for all processes to finish
    for one_process in processes:
        one_process.join()

    # When the processes have finished, combine the outputs together.
    start_time = time.time()
    combine_outputs(header, options)
    end_time = time.time()
    human_time = str(timedelta(seconds=int(end_time - start_time)))
    print("Took %s to combine all outputs." % (human_time))


def getOptions():
    parser = OptionParser()

    parser.add_option("--sam", "-s", dest="sam",
                      help=("Input SAM file containing transcripts to "
                            "correct. Must contain a header."),
                      metavar="FILE", type="string", default="")
    parser.add_option("--genome", "-g", dest="refGenome",
                      help=("Reference genome fasta file. Should be the "
                            "same one used during mapping to generate the "
                            "provided SAM file."),
                      metavar="FILE", type="string", default="")
    parser.add_option("--threads", "-t", dest="n_threads",
                      help=("Number of threads to run program with."),
                      type="int", default=1)
    parser.add_option("--spliceJns", "-j", dest="sjAnnotFile",
                      help=("Splice junction file obtained by mapping "
                            "Illumina reads to the genome using STAR, or "
                            "alternately, extracted from a GTF using the "
                            "accessory script. More formats may be "
                            "supported in the future."),
                      metavar="FILE", type="string", default=None)
    parser.add_option("--variants", "-v", dest="variantFile",
                      help=("VCF formatted file of variants to avoid "
                            "correcting away in the data (optional)."),
                      metavar="FILE", type="string", default=None)
    parser.add_option("--maxLenIndel", dest="maxLenIndel",
                      help="Maximum size indel to correct (Default: 5 bp)",
                      type="int", default=5)
    parser.add_option("--maxSJOffset", dest="maxSJOffset",
                      help=("Maximum distance from annotated splice junction "
                            "to correct (Default: 5 bp)"),
                      type="int", default=5)
    parser.add_option("--outprefix", "-o", dest="outprefix",
                      help=("Output file prefix. '_clean' plus a file "
                            "extension will be added to the end."),
                      metavar="FILE", type="string", default="TC")
    parser.add_option("--correctMismatches", "-m", dest="correctMismatches",
                      help=("If set to false, TranscriptClean will skip "
                            "mismatch correction. Default: true"),
                      type="string", default="true")
    parser.add_option("--correctIndels", "-i", dest="correctIndels",
                      help="If set to false, TranscriptClean will skip indel \
                      correction. Default: true", type="string",
                      default="true")
    parser.add_option("--correctSJs", dest="correctSJs",
                      help=("If set to false, TranscriptClean will skip "
                            "splice junction correction. Default: true, but you must "
                            "provide a splice junction annotation file in order for "
                            "it to work."), type="string", default="true")
    parser.add_option("--dryRun", dest="dryRun", action='store_true',
                      help=("If this option is set, TranscriptClean will "
                            "read in the sam file and record all insertions, "
                            "deletions, and mismatches, but it will skip "
                            "correction. This mode is useful for checking "
                            "the distribution of transcript errors in the "
                            "data before running correction."))
    parser.add_option("--primaryOnly", dest="primaryOnly", action='store_true',
                      help="If this option is set, TranscriptClean will only \
                      output primary mappings of transcripts (ie it will filter \
                      out unmapped and multimapped lines from the SAM input.",
                      default=False)
    parser.add_option("--canonOnly", dest="canonOnly", action='store_true',
                      help=("If this option is set, TranscriptClean will "
                            "output only canonical transcripts and transcripts "
                            "containing annotated noncanonical junctions to the "
                            "clean SAM file at the end of the run."), default=False)
    parser.add_option("--tmpDir", dest="tmp_path",
                      help=("If you would like the tmp files to be written "
                            "somewhere different than the final output, "
                            "provide the path to that location here."),
                      default=None)
    parser.add_option("--bufferSize", dest="buffer_size",
                      help=("Number of lines to output to file at once by "
                            "each thread during run. Default = 100"),
                      type="int", default=100)
    parser.add_option("--deleteTmp", dest="delete_tmp", action='store_true',
                      help=("If this option is set, the temporary directory "
                            "generated by TranscriptClean (TC_tmp) will be "
                            "removed at the end of the run."))

    (options, args) = parser.parse_args()
    options = cleanup_options(options)
    return options


def cleanup_options(options):
    """ Clean up input options by casting to appropriate types etc. """
    options.n_threads = int(options.n_threads)
    options.buffer_size = int(options.buffer_size)
    options.maxLenIndel = int(options.maxLenIndel)
    options.maxSJOffset = int(options.maxSJOffset)
    options.correctIndels = (options.correctIndels).lower()
    options.correctMismatches = (options.correctMismatches).lower()
    options.correctSJs = (options.correctSJs).lower()

    # Screen options for correctness
    if options.correctMismatches not in ["true", "false"]:
        raise RuntimeError(("Invalid choice for --correctMismatches/-m option. "
                            "Valid choices are 'true' and 'false'"))
    if options.correctIndels not in ["true", "false"]:
        raise RuntimeError(("Invalid choice for --correctIndels/-i option. "
                            "Valid choices are 'true' and 'false'"))
    if options.correctSJs not in ["true", "false"]:
        raise RuntimeError(("Invalid choice for --correctSJs option. "
                            "Valid choices are 'true' and 'false'"))

    # Use custom tmp dir location if provided
    if options.tmp_path != None:
        if os.path.isdir(options.tmp_path):
            options.tmp_dir = "/".join(
                (options.tmp_path).split("/") + ["TC_tmp/"])
        else:
            options.tmp_dir = "/".join((options.tmp_path).split("/")
                                       [0:-1] + ["TC_tmp/"])
    else:
        options.tmp_dir = "/".join((options.outprefix).split("/")
                                   [0:-1] + ["TC_tmp/"])

    # If there is a tmp dir there already, remove it
    if os.path.exists(options.tmp_dir):
        os.system("rm -r %s" % options.tmp_dir)
    os.system("mkdir -p %s" % options.tmp_dir)

    # If the specified outprefix is a directory, add TC default prefix to it
    if os.path.isdir(options.outprefix):
        options.outprefix = "/".join((options.outprefix).split("/") + ["TC"])

    return options


def prep_refs(options, transcripts, sam_header):
    """ Process input files and store them in a reference dict.
        SAM transcripts and header are needed if variants provided because
        the variant file will be filtered to include only those that overlap
        the SAM reads.  """

    tmp_dir = options.tmp_dir

    os.system("mkdir -p %s" % tmp_dir)
    genomeFile = options.refGenome
    variantFile = options.variantFile
    sjFile = options.sjAnnotFile

    # Container for references
    refs = Struct()

    # Read in the reference genome.
    print("Reading genome ..............................")
    refs.genome = Fasta(options.refGenome)
    ref_chroms = sorted((refs.genome.keys()))

    # Create a tmp SAM file for the provided reads
    proc = str(os.getpid())
    tmp_sam, sam_chroms = create_tmp_sam(
        sam_header, transcripts, options.tmp_dir, process=proc)

    # Read in splice junctions
    refs.donors = None
    refs.acceptors = None
    refs.sjAnnot = set()
    if options.correctSJs == "false":
        print("SJ correction turned off in options. Skipping reference SJ parsing.")

    elif sjFile != None:
        print("Processing annotated splice junctions ...")
        refs.donors, refs.acceptors, refs.sjAnnot = processSpliceAnnotation(sjFile,
                                                                            tmp_dir, sam_chroms, process=proc)
    else:
        print("No splice annotation provided. Will skip splice junction correction.")

    # Read in variants
    if variantFile != None:
        print("Processing variant file .................")
        refs.snps, refs.insertions, refs.deletions = processVCF(variantFile,
                                                                options.maxLenIndel,
                                                                tmp_dir,
                                                                tmp_sam,
                                                                process=proc)
    else:
        print("No variant file provided. Transcript correction will not be variant-aware.")
        refs["snps"] = refs["insertions"] = refs["deletions"] = {}

    return refs


def create_tmp_sam(sam_header, transcripts, tmp_dir, process="1"):
    """ Put the transcripts in the provided list into a temporary SAM file,
        preceded by the provided header (list form). The file will be located in
        a 'sams' subdir of the tmp_dir provided. Returns the name of the tmp
        file as well as the chromosomes found within"""

    sam_dir = tmp_dir + "split_uncorr_sams/"
    sam_name = sam_dir + process + ".sam"
    os.system("mkdir -p %s" % (sam_dir))

    chroms = set()

    with open(sam_name, 'w') as f:
        for item in sam_header:
            f.write("%s\n" % item)
        for transcript in transcripts:
            f.write("%s\n" % transcript)
            chroms.add(transcript.split('\t')[2])

    return sam_name, chroms


def setup_outfiles(options, process="1"):
    """ Set up output files. If running in parallel, label with a process ID """

    # Place files in a tmp directory
    tmp_dir = options.tmp_dir
    os.system("mkdir -p " + tmp_dir)

    outfiles = Struct()

    # Open sam, fasta, and log  outfiles
    oSam = open(tmp_dir + "clean_" + process + ".sam", 'w')
    oFa = open(tmp_dir + "clean_" + process + ".fa", 'w')
    transcriptLog = open(tmp_dir + "clean_" + process + ".log", 'w')
    transcriptErrorLog = open(tmp_dir + "clean_" + process + ".TElog", 'w')

    outfiles.sam = oSam
    outfiles.fasta = oFa
    outfiles.log = transcriptLog
    outfiles.TElog = transcriptErrorLog

    return outfiles


def close_outfiles(outfiles):
    """ Close all of the output files """

    for f in list(outfiles.values()):
        f.close()

    return


def transcript_init(transcript_line, genome, sjAnnot):
    """ Attempt to initialize a Transcript object from a SAM entry. First,
        create a log object, and note the mapping status of the read. If the
        transcript alignment is unmapped or non-primary, then create a log
        entry for the transcript, but do not bother initializing a transcript
        object. output that Otherwise, return the transcript object along with
        the log.

        Input:
            - Transcript SAM line (string)
            - Reference genome (needed to determine canonical-ness of jns
            - Splice annot set object (needed to look up whether a junction
              is annotated or not)

        Returns:
            - Transcript object if primary alignment, None if not
            - logInfo object
     """

    # Check mapping
    sam_fields = transcript_line.split('\t')
    logInfo = init_log_info(sam_fields)

    if logInfo.Mapping != "primary":
        return None, logInfo

    # try:
    transcript = Transcript(sam_fields, genome, sjAnnot)

    # except Exception as e:
    #     warnings.warn("Problem parsing transcript with ID '" +
    #                   logInfo.TranscriptID + "'")
    #     print(e)
    #     return None, logInfo

    return transcript, logInfo


def batch_correct(sam_transcripts, options, refs, outfiles, buffer_size=100):
    """Correct and output n lines of transcript SAM file at a time in a batch.
       The purpose is to stagger when the different threads are writing to disk.
       For a given transcript, there can be only one sam, fasta, and log entry,
       but there can be more than one transcript error (TE)."""

    print("Correcting transcripts...")

    outSam = outfiles.sam
    outFa = outfiles.fasta
    tL = outfiles.log
    tE = outfiles.TElog
    sam_lines = fasta_lines = log_lines = TE_log_lines = ""

    counter = 0
    for transcript_line in sam_transcripts:
        transcript, logInfo, TE_entries = correct_transcript(transcript_line,
                                                             options, refs)

        # If the transcript object returned is None, this means that the
        # current read was not a primary mapper. Add to the log file,
        # but do not output a fasta sequence. Only output the original
        # SAM alignment if primaryOnly and canonOnly modes are off
        if transcript == None:
            log_lines += create_log_string(logInfo) + "\n"
            if options.primaryOnly == False and options.canonOnly == False:
                sam_lines += transcript_line + "\n"

        # Output the transcript in SAM and Fasta format, plus output the logs
        else:
            canonical_or_annot = transcript.isCanonical or transcript.allJnsAnnotated
            if options.canonOnly == False or canonical_or_annot == True:
                sam_lines += transcript.printableSAM() + "\n"
                fasta_lines += transcript.printableFa() + "\n"
            log_lines += create_log_string(logInfo) + "\n"
            TE_log_lines += TE_entries

        # Check whether to empty the buffer
        if counter > buffer_size:
            outSam.write(sam_lines)
            outFa.write(fasta_lines)
            tL.write(log_lines)
            tE.write(TE_log_lines)
            sam_lines = fasta_lines = log_lines = TE_log_lines = ""
            counter = 0
        else:
            counter += 1

    # Write any remaining lines
    outSam.write(sam_lines)
    outFa.write(fasta_lines)
    tL.write(log_lines)
    tE.write(TE_log_lines)
    return


def correct_transcript(transcript_line, options, refs):
    """ Given a line from a SAM file, create a transcript object. If it's a
        primary alignment, then perform the corrections specified in the
        options.
    """
    orig_transcript, logInfo = transcript_init(transcript_line, refs.genome,
                                               refs.sjAnnot)
    TE_entries = ""

    if orig_transcript == None:
        return orig_transcript, logInfo, TE_entries

    # Correct the transcript
    try:
        upd_transcript = copy(orig_transcript)
        upd_logInfo = logInfo
        # Mismatch correction
        if options.correctMismatches == "true":
            mismatch_TE = correctMismatches(upd_transcript, refs.genome,
                                            refs.snps, upd_logInfo)
            if mismatch_TE != "":
                TE_entries += mismatch_TE

        if options.correctIndels == "true":
            # Insertion correction
            ins_TE = correctInsertions(upd_transcript, refs.genome, refs.insertions,
                                       options.maxLenIndel, upd_logInfo)
            if ins_TE != "":
                TE_entries += ins_TE

            # Deletion correction
            del_TE = correctDeletions(upd_transcript, refs.genome, refs.deletions,
                                      options.maxLenIndel, upd_logInfo)
            if del_TE != "":
                TE_entries += del_TE

        # NCSJ correction
        if len(refs.sjAnnot) > 0 and options.correctSJs == "true":
            upd_transcript, ncsj_TE = cleanNoncanonical(upd_transcript, refs,
                                                        options.maxSJOffset,
                                                        upd_logInfo)
            if ncsj_TE != "":
                TE_entries += ncsj_TE

    except Exception as e:
        warnings.warn(("Problem encountered while correcting transcript "
                       "with ID %s. Will output original version.") %
                      orig_transcript.QNAME)
        print(e)
        TE_entries = ""
        return orig_transcript, logInfo, TE_entries

    # After successful correction, return transcript object, updated
    # logInfo, and the TE log lines generated during correction
    return upd_transcript, upd_logInfo, TE_entries


def validate_chroms(genome_file, vcf_file, sam_chroms):
    """ Make sure that every chromosome in the SAM file also exists in the
        reference genome. This is a common source of crashes. Also, make sure
        that the VCF chromosome naming convention matches the SAM file if the
        user has provided a VCF. """

    genome = Fasta(genome_file)
    fasta_chroms = set(genome.keys())

    # Remove '*' from sam chromosome set if present
    if "*" in sam_chroms:
        sam_chroms.remove("*")

    # Check whether all of the sam chromosomes are in the fasta file.
    # If not, raise an error
    if not sam_chroms.issubset(fasta_chroms):
        sam_chroms = "{" + \
            ", ".join(['"' + str(x) + '"' for x in sam_chroms]) + '}'
        fasta_chroms = "{" + \
            ", ".join(['"' + str(x) + '"' for x in fasta_chroms]) + '}'
        error_msg = 'One or more SAM chromosomes were not found in the fasta reference.\n' + \
                    'SAM chromosomes:\n' + sam_chroms + '\n' + \
                    'FASTA chromosomes:\n' + fasta_chroms + '\n' + \
                    ("One common cause of this problem is when the fasta headers "
                     "contain more than one word. If this is the case, try "
                     "trimming the headers to include only the chromosome name "
                     "(i.e. '>chr1').")
        raise RuntimeError(error_msg)

    # Check VCF chroms
    if vcf_file != None:
        try:
            import pybedtools
            vcf_obj = pybedtools.BedTool(vcf_file)
        except Exception as e:
            print(e)
            raise RuntimeError(("Problem reading the provided VCF file. Please "
                                "make sure that the filename is correct, that "
                                "pybedtools is installed, and that your VCF "
                                "file is properly formatted, including a header."))

        vcf_chroms = set([x.chrom for x in vcf_obj])
        if len(sam_chroms.intersection(vcf_chroms)) == 0:
            sam_chroms = "{" + \
                ", ".join(['"' + str(x) + '"' for x in sam_chroms]) + '}'
            vcf_chroms = "{" + \
                ", ".join(['"' + str(x) + '"' for x in vcf_chroms]) + '}'
            error_msg = ("None of the chromosomes included in the VCF file "
                         "matched the chromosomes in the provided SAM file.\n"
                         'SAM chromosomes:\n' + sam_chroms + '\n'
                         'VCF chromosomes:\n' + vcf_chroms + '\n')
            raise RuntimeError(error_msg)

    return


def split_SAM(samFile, n):
    """ Given a sam file, separate the header line from the remaining entries,
        record which chromosomes are present in the reads, and split the reads
        into n chunks.
    """
    header = []
    chroms = set()
    transcript_lines = []

    with open(samFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('@'):
                header.append(line)
            else:
                transcript_lines.append(line)
                chrom = line.split("\t")[2]
                chroms.add(chrom)

    # Now split the transcripts into n chunks
    chunks = split_input(transcript_lines, n)

    return header, chroms, chunks


def split_input(my_list, n):
    """ Splits input into n sublists of roughly equal size"""

    chunks = []
    index = 0
    batch_size = ceil(len(my_list)/n)
    while index < len(my_list):
        try:
            batch = my_list[index:index + batch_size]
        except:
            batch = my_list[index:]
        chunks.append(batch)
        index += batch_size
    return chunks


def run_chunk(transcripts, options, sam_header):
    """ Contains everything needed to run a subset of the transcript input
        on its own core """

    # Prep the references (i.e. genome, sjs, variants)
    start_time = time.time()
    refs = prep_refs(options, transcripts, sam_header)
    end_time = time.time()
    time_slice = end_time - start_time
    human_time = str(timedelta(seconds=int(time_slice)))
    print("Reference file processing took %s" % (human_time))

    # Set up the outfiles
    outfiles = setup_outfiles(options, str(os.getpid()))

    # Correct the transcripts
    start_time = time.time()
    batch_correct(transcripts, options, refs, outfiles,
                  buffer_size=options.buffer_size)
    end_time = time.time()
    human_time = str(timedelta(seconds=int(end_time - start_time)))
    print("Took %s to process transcript batch." % (human_time))

    # Close up the outfiles
    close_outfiles(outfiles)

    return


def run_chunk_dryRun(transcripts, options):
    """ Run dryRun mode on a set of transcripts on one core """

    # Set up the outfiles
    outfiles = setup_outfiles(options, str(os.getpid()))

    # Run dryRun mode
    dryRun(transcripts, options, outfiles)

    # Close up the outfiles
    close_outfiles(outfiles)

    return


def combine_outputs(sam_header, options):
    """ Combine sam, fasta, and log outputs from separate sub-runs into the
        final output files. Then, remove the temporary directory.
    """

    outprefix = options.outprefix
    tmp_dir = options.tmp_dir

    # Make filenames
    sam = outprefix + "_clean.sam"
    fasta = outprefix + "_clean.fa"
    log = outprefix + "_clean.log"
    TElog = outprefix + "_clean.TE.log"

    # Open log outfiles
    if os.path.exists(log):
        os.remove(log)
    if os.path.exists(TElog):
        os.remove(TElog)

    transcriptLog = open(log, 'w')
    transcriptErrorLog = open(TElog, 'w')

    # Add headers to logs
    transcriptLog.write("\t".join(["TranscriptID", "Mapping",
                                   "corrected_deletions", "uncorrected_deletions",
                                   "variant_deletions", "corrected_insertions",
                                   "uncorrected_insertions", "variant_insertions",
                                   "corrected_mismatches", "uncorrected_mismatches",
                                   "corrected_NC_SJs", "uncorrected_NC_SJs"]) + "\n")

    transcriptErrorLog.write("\t".join(["TranscriptID", "Position", "ErrorType",
                                        "Size", "Corrected", "ReasonNotCorrected"]) + "\n")

    # Add header to sam file
    if options.dryRun != True:
        if os.path.exists(sam):
            os.remove(sam)

        oSam = open(sam, 'w')
        for line in sam_header:
            oSam.write(line + "\n")
        oSam.close()

    # Close files
    transcriptLog.close()
    transcriptErrorLog.close()

    # Now combine the subfiles
    if options.dryRun != True:
        os.system("cat %s/*.sam >> %s" % (tmp_dir, sam))
        os.system("cat %s/*.fa > %s" % (tmp_dir, fasta))
    os.system("cat %s/*.log >> %s" % (tmp_dir, log))
    os.system("cat %s/*.TElog >> %s" % (tmp_dir, TElog))

    # Clean up the temporary directory if requested
    if options.delete_tmp:
        os.system("rm -r %s" % tmp_dir)

    return


def processSpliceAnnotation(annotFile, tmp_dir, read_chroms, process="1"):
    """ Reads in the tab-separated STAR splice junction file and creates a
        bedtools object. Also creates a dict (annot) to allow easy lookup
        to find out if a splice junction is annotated or not. Only junctions
        located on the provided chromosomes are included."""

    bedstr = ""
    annot = set()
    tmp_dir = tmp_dir + "splice_files/"
    os.system("mkdir -p %s" % (tmp_dir))

    donor_file = tmp_dir + "%s_ref_splice_donors_tmp.bed" % (process)
    acceptor_file = tmp_dir + "%s_ref_splice_acceptors_tmp.bed" % (process)
    o_donor = open(donor_file, 'w')
    o_acceptor = open(acceptor_file, 'w')
    with open(annotFile, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            if chrom not in read_chroms:
                continue

            start = int(fields[1])
            end = int(fields[2])

            strand = "."
            if fields[3] == "1":
                strand = "+"
            if fields[3] == "2":
                strand = "-"

            intronMotif = int(fields[4])
            annotated = int(fields[5])
            uniqueReads = fields[6]
            multiReads = fields[7]
            maxOverhang = fields[8]

            # ID the splice donor/acceptor identity
            if strand == "+":
                file1 = o_donor
                file2 = o_acceptor
                type1 = "donor"
                type2 = "acceptor"
            elif strand == "-":
                file1 = o_acceptor
                file2 = o_donor
                type1 = "acceptor"
                type2 = "donor"
            else:
                continue

            # Make one bed entry for each end of the junction and write to
            # splice donor and acceptor files
            bed1 = "\t".join(
                [chrom, str(start - 1), str(start), ".", "name", strand])
            bed2 = "\t".join(
                [chrom, str(end - 1), str(end), ".", "name", strand])
            file1.write(bed1 + "\n")
            file2.write(bed2 + "\n")

            # Add a record of the entire junction to the annot set to indicate
            # that it is annotated
            junction = "_".join([chrom, str(start), strand]) + "," + \
                       "_".join([chrom, str(end), strand])
            annot.add(junction)
    o_donor.close()
    o_acceptor.close()

    # Convert bed files into BedTool objects
    # donor_sorted = tmp_dir + "%s_ref_splice_donors_tmp.sorted.bed" % (process)
    # acceptor_sorted = tmp_dir + \
    #     "%s_ref_splice_acceptors_tmp.sorted.bed" % (process)
    # os.system('sort -u %s | bedtools sort -i - > %s' %
    #           (donor_file, donor_sorted))
    # os.system('sort -u %s | bedtools sort -i - > %s' %
    #           (acceptor_file, acceptor_sorted))

    # splice_donor_bedtool = pybedtools.BedTool(donor_sorted)
    # splice_acceptor_bedtool = pybedtools.BedTool(acceptor_sorted)

    import pyranges as pr
    # import pdb
    # pdb.set_trace()

    try:
        splice_donor = pr.read_bed(donor_file).drop_duplicate_positions()
    except IndexError:
        splice_donor = pr.PyRanges()

    try:
        splice_acceptor = pr.read_bed(acceptor_file).drop_duplicate_positions()
    except IndexError:
        splice_acceptor = pr.PyRanges()

    # Raise a warning if there are no splice donors or acceptors
    if splice_donor.length == 0 or splice_acceptor.length == 0:
        warnings.warn(("Warning: No splice donors or acceptors found on "
                       "chromosomes: %s. If this is unexpected, check your SJ "
                       "annotation file." % (", ".join(read_chroms))))

    return splice_donor, splice_acceptor, annot


def processVCF(vcf, maxLen, tmp_dir, sam_file, process="1"):
    """ This function reads in variants from a VCF file and stores them. SNPs
        are stored in a dictionary by position. Indels are stored in their own
        dictionary by start and end position. A pre-filtering step ensures that
        variants are included only if they overlap with the input SAM reads.
        This is a space-saving measure. If add_chr is set to 'True', the prefix
        'chr' will be added to each chromosome name."""

    SNPs = {}
    insertions = {}
    deletions = {}

    # Create a temporary BAM version of the input reads
    tmp_bam_dir = tmp_dir + "split_uncorr_bams/"
    os.system("mkdir -p %s" % (tmp_bam_dir))
    bam = tmp_bam_dir + "%s_tmp_reads.bam" % (process)
    try:
        os.system("samtools view -bS %s > %s" % (sam_file, bam))

    except Exception as e:
        print(e)
        raise RuntimeError(("Problem converting SAM file to BAM for variant "
                            "prefiltering step. Please make sure that you have "
                            "Samtools installed, and that your SAM "
                            "file is properly formatted, including a header."))

    # Create a BedTool object for the VCF
    try:
        import pybedtools
        vcf_obj = pybedtools.BedTool(vcf)
    except Exception as e:
        print(e)
        raise RuntimeError(("Problem reading the provided VCF file. Please make "
                            "sure that Bedtools is installed, and that your VCF "
                            "file is properly formatted, including a header. "))

    # Filter the VCF file to include only those variants that overlap the reads
    filtered_vcf = vcf_obj.intersect(bam, u=True, nonamecheck=True)

    if len(filtered_vcf) == 0:
        warnings.warn(("Warning: none of variants provided overlapped the"
                       "input reads."))
        return SNPs, insertions, deletions

    # Now parse the filtered VCF
    for variant in filtered_vcf:

        chrom = variant[0]
        pos = int(variant[1])
        ref = variant[3]
        alt = variant[4].split(",")
        refLen = len(ref)

        # It is possible for a single position to have more than one
        # type of alternate allele (ie mismatch or indel). So it is
        # necessary to treat each alternate allele separately.
        for allele in alt:
            altLen = len(allele)

            # SNP case
            if refLen == altLen == 1:
                ID = chrom + "_" + str(pos)
                if ID not in SNPs:
                    SNPs[ID] = [allele]
                else:
                    SNPs[ID].append(allele)
            # Insertion/Deletion
            else:
                size = abs(refLen - altLen)
                # Only store indels of correctable size
                if size > maxLen:
                    continue
                # Positions in VCF files are one-based.
                # Inserton/deletion sequences always start with the
                # preceding normal reference base
                actPos = pos + 1
                allele = allele[1:]
                ID = "_".join([chrom, str(pos), str(actPos + size - 1)])

                if refLen - altLen < 0:  # Insertion
                    if ID not in insertions:
                        insertions[ID] = [allele]
                    else:
                        insertions[ID].append(allele)
                elif refLen - altLen > 0:  # Deletion
                    deletions[ID] = 1

    return SNPs, insertions, deletions


def correctInsertions(transcript, genome, variants, maxLen, logInfo):
    """ Corrects insertions up to size maxLen using the reference genome.
        If a variant file was provided, correction will be SNP-aware. """

    logInfo.uncorrected_insertions = 0
    logInfo.corrected_insertions = 0
    logInfo.variant_insertions = 0
    TE_entries = ""

    origSeq = transcript.SEQ
    origCIGAR = transcript.CIGAR
    transcript_ID = transcript.QNAME

    cigarOps, cigarCounts = transcript.splitCIGAR()

    # Check for insertions. If none are present, we can skip this transcript
    if "I" not in origCIGAR:
        return TE_entries

    newCIGAR = ""
    newSeq = ""
    MVal = 0
    seqPos = 0

    # Start at position in the genome where the transcript starts.
    genomePos = transcript.POS

    # Iterate over operations to sequence and repair insertions
    for op, ct in zip(cigarOps, cigarCounts):

        currPos = transcript.CHROM + ":" + str(genomePos) + "-" + \
            str(genomePos + ct - 1)

        if op == "M":
            newSeq = newSeq + origSeq[seqPos:seqPos + ct]
            MVal += ct
            seqPos += ct
            genomePos += ct

        if op == "I":
            ID = "_".join([transcript.CHROM, str(genomePos - 1),
                           str(genomePos + ct - 1)])

            # Only insertions of a given size are corrected
            if ct <= maxLen:
                # Check if the insertion is in the optional variant catalog.
                if ID in variants:
                    # The insertion perfectly matches a variant position.
                    # Leave the sequence alone if it matches an allele sequence.
                    currSeq = origSeq[seqPos:seqPos + ct]
                    if currSeq in variants[ID]:
                        logInfo.variant_insertions += 1
                        errorEntry = "\t".join([transcript_ID, ID, "Insertion",
                                                str(ct), "Uncorrected",
                                                "VariantMatch"])
                        TE_entries += errorEntry + "\n"

                        # Leave insertion in
                        MVal, newCIGAR = endMatch(MVal, newCIGAR)
                        newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                        newCIGAR = newCIGAR + str(ct) + op
                        seqPos += ct
                        continue

                # Correct insertion
                errorEntry = "\t".join([transcript_ID, ID, "Insertion", str(ct),
                                        "Corrected", "NA"])
                logInfo.corrected_insertions += 1
                TE_entries += errorEntry + "\n"

                # Subtract the inserted bases by skipping them.
                # GenomePos stays the same, as does MVal
                seqPos += ct
            else:  # Move on without correcting insertion because it is too big
                errorEntry = "\t".join([transcript_ID, ID, "Insertion", str(ct),
                                        "Uncorrected", "TooLarge"])
                logInfo.uncorrected_insertions += 1
                TE_entries += errorEntry + "\n"
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

        # N, H, and D operations are cases where the transcript sequence
        # is missing bases that are in the reference genome.
        if op in ["N", "H", "D"]:
            # End any ongoing match
            MVal, newCIGAR = endMatch(MVal, newCIGAR)
            genomePos += ct
            newCIGAR = newCIGAR + str(ct) + op

    # End any ongoing match
    MVal, newCIGAR = endMatch(MVal, newCIGAR)

    # Update transcript
    transcript.CIGAR = newCIGAR
    transcript.SEQ = newSeq

    return TE_entries


def correctDeletions(transcript, genome, variants, maxLen, logInfo):
    """ Corrects deletions up to size maxLen using the reference genome.
        If a variant file was provided, correction will be variant-aware."""

    logInfo.uncorrected_deletions = 0
    logInfo.corrected_deletions = 0
    logInfo.variant_deletions = 0
    TE_entries = ""

    transcript_ID = transcript.QNAME
    chrom = transcript.CHROM
    origSeq = transcript.SEQ
    origCIGAR = transcript.CIGAR

    cigarOps, cigarCounts = transcript.splitCIGAR()

    # Check for deletions. If none are present, we can skip this transcript
    if "D" not in origCIGAR:
        return TE_entries

    newCIGAR = ""
    newSeq = ""
    MVal = 0
    seqPos = 0

    # Start at position in the genome where the transcript starts.
    genomePos = transcript.POS

    # Iterate over operations to sequence and repair mismatches and microindels
    for op, ct in zip(cigarOps, cigarCounts):

        currPos = chrom + ":" + str(genomePos) + "-" + \
            str(genomePos + ct - 1)

        if op == "M":
            newSeq = newSeq + origSeq[seqPos:seqPos + ct]
            MVal += ct
            seqPos += ct
            genomePos += ct

        if op == "D":
            ID = "_".join([chrom, str(genomePos - 1), str(genomePos + ct - 1)])
            if ct <= maxLen:

                # Check if the deletion is in the optional variant catalog.
                if ID in variants:
                    # The deletion perfectly matches a deletion variant. Leave the deletion in.
                    # currSeq = genome.sequence(
                    #     {'chr': chrom, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                    currSeq = genome.get_seq(chrom, genomePos, genomePos+ct-1).seq
                    errorEntry = "\t".join(
                        [transcript_ID, ID, "Deletion", str(ct), "Uncorrected", "VariantMatch"])
                    logInfo.variant_deletions += 1
                    TE_entries += errorEntry + "\n"

                    MVal, newCIGAR = endMatch(MVal, newCIGAR)
                    genomePos += ct
                    newCIGAR = newCIGAR + str(ct) + op
                    continue

                # Correct deletion if we're not in variant-aware mode
                errorEntry = "\t".join(
                    [transcript_ID, ID, "Deletion", str(ct), "Corrected", "NA"])
                logInfo.corrected_deletions += 1
                TE_entries += errorEntry + "\n"

                # Add the missing reference bases
                # refBases = genome.sequence(
                #     {'chr': chrom, 'start': genomePos, 'stop': genomePos + ct - 1}, one_based=True)
                refBases = genome.get_seq(chrom, genomePos, genomePos+ct-1).seq
                newSeq = newSeq + refBases
                genomePos += ct
                MVal += ct

            # Deletion is too big to fix
            else:
                errorEntry = "\t".join(
                    [transcript_ID, ID, "Deletion", str(ct), "Uncorrected", "TooLarge"])
                logInfo.uncorrected_deletions += 1
                TE_entries += errorEntry + "\n"

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

    transcript.CIGAR = newCIGAR
    transcript.SEQ = newSeq
    transcript.NM, transcript.MD = transcript.getNMandMDFlags(genome)
    return TE_entries


def correctMismatches(transcript, genome, variants, logInfo):
    """ This function corrects mismatches in the provided transcript. If a
        variant file was provided, correction will be variant-aware."""

    logInfo.uncorrected_mismatches = 0
    logInfo.corrected_mismatches = 0
    TE_entries = ""

    origSeq = transcript.SEQ
    origCIGAR = transcript.CIGAR
    origMD = transcript.MD

    # Check for mismatches. If none are present, we can skip this transcript
    if any(i in origMD.upper() for i in 'ACTGN') == False:
        return TE_entries

    newCIGAR = ""
    newSeq = ""
    MVal = 0
    seqPos = 0
    genomePos = transcript.POS

    # Merge CIGAR and MD tag information so that we have the locations of all
    # insertions, deletions, and mismatches
    mergeOperations, mergeCounts = transcript.mergeMDwithCIGAR()

    # Iterate over operations to sequence and repair mismatches and microindels
    for op, ct in zip(mergeOperations, mergeCounts):
        if op == "M":
            newSeq = newSeq + origSeq[seqPos:seqPos + ct]
            MVal += ct
            seqPos += ct
            genomePos += ct

        # This denotes a mismatch
        if op == "X":
            # Check if the position matches a variant in the optional variant catalog.
            ID = transcript.CHROM + "_" + str(genomePos)
            if ID in variants:
                # Mismatch position matches a known variant. If the base matches an alternative allele, do not correct it.
                currBase = origSeq[seqPos]
                if currBase in variants[ID]:
                    errorEntry = "\t".join([transcript.QNAME, ID, "Mismatch",
                                            str(ct), "Uncorrected",
                                            "VariantMatch"])
                    TE_entries += errorEntry + "\n"
                    logInfo.uncorrected_mismatches += 1

                    # Keep the base as-is
                    newSeq = newSeq + origSeq[seqPos:seqPos + ct]
                    MVal += ct
                    seqPos += ct
                    genomePos += ct
                    continue
            # Otherwise, correct the mismatch to reference base
            errorEntry = "\t".join([transcript.QNAME, ID, "Mismatch",
                                    str(ct), "Corrected", "NA"])
            TE_entries += errorEntry + "\n"
            logInfo.corrected_mismatches += 1

            # Change sequence base to the reference base at this position
            # newSeq = newSeq + genome.sequence({'chr': transcript.CHROM,
            #                                    'start': genomePos,
            #                                    'stop': genomePos + ct - 1},
            #                                   one_based=True)
            newSeq += genome.get_seq(transcript.CHROM, genomePos, genomePos+ct-1).seq
            seqPos += ct  # skip the original sequence base
            genomePos += ct  # advance the genome position
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

    transcript.CIGAR = newCIGAR
    transcript.SEQ = newSeq
    # TODO: find faster way to compute
    transcript.NM, transcript.MD = transcript.getNMandMDFlags(genome)

    return TE_entries


def endMatch(MVal, newCIGAR):
    """ Function to end an ongoing match during error correction, updating new
        CIGAR and MD tag strings as needed. MVal keeps track of how many
        matched bases we have at the moment so that when errors are fixed,
        adjoining matches can be combined. When an intron or unfixable error
        is encountered, it is necessary to end the match by outputting the
        current MVal and then setting the variable to zero."""

    if MVal > 0:
        newCIGAR = newCIGAR + str(MVal) + "M"
        MVal = 0

    return MVal, newCIGAR


def cleanNoncanonical(transcript, refs, maxDist, logInfo):
    """ Iterate over each junction in the provided transcript object. If the
        junction is noncanonical, attempt to correct it, and record the result.
        Returns a modified transcript object that originates from a copy of the
        object passed into the function, along with the error log entries.
    """
    logInfo.corrected_NC_SJs = 0
    logInfo.uncorrected_NC_SJs = 0
    TE_entries = ""

    if transcript.isCanonical == True:
        return transcript, TE_entries

    # Iterate over all junctions and attempt to correct
    for splice_jn_num in range(0, len(transcript.spliceJunctions)):
        junction = transcript.spliceJunctions[splice_jn_num]
        if junction.isCanonical:
            continue
        else:
            temp_transcript = copy(transcript)
            ID = "_".join([junction.chrom, str(
                junction.start), str(junction.end)])
            correction_status, reason, dist = attempt_jn_correction(temp_transcript,
                                                                    splice_jn_num,
                                                                    refs.genome,
                                                                    refs.donors,
                                                                    refs.acceptors,
                                                                    refs.sjAnnot,
                                                                    maxDist)

            # If there were no problems during correction, replace the original
            # transcript with the modified copy.
            if correction_status == True:
                transcript = temp_transcript
                logInfo.corrected_NC_SJs += 1
                status = "Corrected"
            else:
                logInfo.uncorrected_NC_SJs += 1
                status = "Uncorrected"

            # Update transcript error log
            errorEntry = "\t".join([transcript.QNAME, ID, "NC_SJ_boundary",
                                    str(dist), status, reason])
            TE_entries += errorEntry + "\n"

    return transcript, TE_entries


def find_closest_bound(sj_bound, ref_bounds):
    """ Given one side of a splice junction, find the closest reference """

    import pyranges as pr

    # gr = pr.read_bed(ref_bounds.fn)

    row = pr.PyRanges(chromosomes=[sj_bound.chrom],
                      starts=[sj_bound.pos - 1],
                      ends=[sj_bound.pos],
                      strands=[sj_bound.strand])
    match = row.nearest(ref_bounds).df.iloc[0]
    # Create a Bedtool object for the bound
    # bed_pos = pybedtools.BedTool(sj_bound.getBED(), from_string=True)

    # # Run Bedtools Closest operation
    # closest = bed_pos.closest(
    #     ref_bounds, s=True, D="ref", t="first", nonamecheck=True)[0]

    # Create an object to represent the closest match
    # Coordinates are 0-based since they are coming from BED
    obj_closest = Struct()
    obj_closest.chrom = match.Chromosome  # closest[6]
    obj_closest.start = match.Start_b  # int(closest[7])
    obj_closest.end = match.End_b  # int(closest[8])
    obj_closest.dist = match.Start_b - match.Start  # int(closest[-1])

    return obj_closest


def find_closest_ref_junction(junction, ref_donors, ref_acceptors):
    """ Given a splice junction object and a splice reference, locate the
        closest junction in the provided splice reference BedTool object"""

    # Get the splice donor and acceptor of the junction
    donor = junction.get_splice_donor()
    acceptor = junction.get_splice_acceptor()

    # Find the closest match for each
    closest_ref_donor = find_closest_bound(donor, ref_donors)
    closest_ref_acceptor = find_closest_bound(acceptor, ref_acceptors)

    return closest_ref_donor, closest_ref_acceptor


def update_post_ncsj_correction(transcript, splice_jn_num, genome, sjAnnot):
    """ After correcting a noncanonical splice junction, perform the
        following updates:
        - Reassign the position of the junction (based on its bounds)
        - Recompute the splice motif and assign to junction
        - Recompute NM/MD tags
        - Check canonical and annotated status
    """
    junction = transcript.spliceJunctions[splice_jn_num]
    junction.recheckPosition()
    junction.checkSpliceMotif(genome, sjAnnot)
    transcript.NM, transcript.MD = transcript.getNMandMDFlags(genome)
    transcript.jM, transcript.jI = transcript.get_jM_jI_tags_from_sjs()
    transcript.isCanonical = transcript.recheckCanonical()
    transcript.allJnsAnnotated = transcript.recheckJnsAnnotated()
    return


def attempt_jn_correction(transcript, splice_jn_num, genome, ref_donors,
                          ref_acceptors, sjAnnot, maxDist):
    """ Given a noncanonical splice junction, try to correct it by locating
        a nearby annotated junction within the allowable distance.

        Returns:
            Corrected (True/False)
            Reason for no correction ("NA"/"TooFarFromAnnotJn"/"Other")
            CombinedDist
     """

    junction = transcript.spliceJunctions[splice_jn_num]

    # Find the closest reference junction
    ref_donor, ref_acceptor = find_closest_ref_junction(junction, ref_donors,
                                                        ref_acceptors)

    # Compute the overall distance of the original junction from the reference
    combined_dist = combinedJunctionDist(ref_donor.dist, ref_acceptor.dist)

    # Only attempt to rescue junction boundaries that are within maxDist bp of
    # an annotated junction
    if combined_dist > maxDist or abs(ref_donor.dist) > 2*maxDist or \
            abs(ref_acceptor.dist) > 2*maxDist:
        return False, "TooFarFromAnnotJn", combined_dist

    try:
        # Attempt to fix the splice donor side
        donor = junction.get_splice_donor()
        transcript.SEQ, transcript.CIGAR = fix_one_side_of_junction(transcript.CHROM,
                                                                    transcript.POS, splice_jn_num,
                                                                    donor, ref_donor.dist, genome,
                                                                    transcript.SEQ, transcript.CIGAR)

        # Attempt to fix the splice acceptor side
        acceptor = junction.get_splice_acceptor()
        transcript.SEQ, transcript.CIGAR = fix_one_side_of_junction(transcript.CHROM,
                                                                    transcript.POS, splice_jn_num,
                                                                    acceptor, ref_acceptor.dist, genome,
                                                                    transcript.SEQ, transcript.CIGAR)

        # Now, perform updates:
        update_post_ncsj_correction(transcript, splice_jn_num, genome, sjAnnot)

    except Exception as e:
        return False, "Other", combined_dist

    return True, "NA", combined_dist


def combinedJunctionDist(dist_0, dist_1):
    """Computes the combined genomic distance of two splice junction ends
    from the closest annotated junctions. In essence, it is finding the
    size indel that could have created the discrepancy between the
    reference and transcript junctions.

    Examples ('|' character respresents end of exon):

        Reference:     ----->|          |<-----
        Transcript:      ----->|      |<-----
            dist_0 = -2, dist_1 = +2, combined dist = 4

        Reference:     ----->|          |<-----
        Transcript:    ----->|        |<-----
            dist_0 = 0, dist_1 = +2, combined dist = 2

        Reference:     ----->|          |<-----
        Transcript:   ----->|       |<-----
            dist_0 = +1, dist_1 = +4, combined dist = 3

    """
    # If dist_0 and dist_1 have different signs, the combined distance is
    # the sum of their absolute values
    # https://stackoverflow.com/questions/66882/simplest-way-to-check-if-two-integers-have-same-sign
    if ((dist_0<0)!=(dist_1<0)):
    # if dist_0*dist_1 <= 0:
        combined_dist = abs(dist_0) + abs(dist_1)
    else:
        combined_dist = abs(abs(dist_0) - abs(dist_1))
    return combined_dist


def group_CIGAR_by_exon_and_intron(CIGAR, seq):
    """ Given a CIGAR string, group all of the operations by exon or intron,
        resulting in one string per exon/intron. Also, subdivide the provided
        sequence by exon as well."""

    # Inits
    currExonStr = ""
    currExonCIGAR = ""
    currSeqIndex = 0
    operations, counts = splitCIGAR(CIGAR)

    # Outputs
    exonSeqs = []
    exonCIGARs = []
    intronCIGARs = []

    for op, ct in zip(operations, counts):
        if op == "N":
            exonSeqs.append(currExonStr)
            intronCIGARs.append(ct)
            exonCIGARs.append(currExonCIGAR)
            currExonStr = ""
            currExonCIGAR = ""
        elif op in ["M", "S", "I"]:
            currExonStr = currExonStr + str(seq[currSeqIndex:currSeqIndex+ct])
            currExonCIGAR = currExonCIGAR + str(ct) + op
            currSeqIndex += ct
        else:
            currExonCIGAR = currExonCIGAR + str(ct) + op
    # Append last exon entries
    exonSeqs.append(currExonStr)
    exonCIGARs.append(currExonCIGAR)

    return exonCIGARs, intronCIGARs, exonSeqs


def fix_one_side_of_junction(chrom, transcript_start, jn_number, intronBound, d,
                             genome, seq, oldCIGAR):
    """ Corrects one side of a noncanonical splice junction by the specified
        number of bases. This involves modifying the sequence and CIGAR string.

        Returns modified CIGAR string and sequence
    """

    # If d = 0, then the bound already matches the annotation.
    if d == 0:
        return seq, oldCIGAR

    # First, split CIGAR by exon/intron, as well as the sequence
    exonCIGARs, intronCIGARs, exonSeqs = group_CIGAR_by_exon_and_intron(
        oldCIGAR, seq)

    # Now go and fix the specific splice junction piece
    if intronBound.bound == 0:
        targetExon = jn_number
        exon = exonSeqs[targetExon]
        if d > 0:  # Need to add d bases from reference to end of exon. Case 1.
            # For CIGAR string,
            exonEnd = intronBound.pos - 1
            seqIndex = exonEnd - transcript_start + 1
            # refAdd = genome.sequence(
            #     {'chr': chrom, 'start': exonEnd + 1, 'stop': exonEnd + d}, one_based=True)
            refAdd = genome.get_seq(chrom, exonEnd+1, exonEnd+d).seq
            exonSeqs[targetExon] = exon + refAdd
            intronBound.pos += d
        if d < 0:  # Need to subtract from end of exon sequence. Case 3
            exonSeqs[targetExon] = exon[0:d]
            intronBound.pos += d
        intronCIGARs[jn_number] -= d
        exonCIGARs[targetExon] = editExonCIGAR(exonCIGARs[targetExon], -1, d)
    else:
        targetExon = jn_number + 1
        exon = exonSeqs[targetExon]
        if d < 0:  # Need to add d bases from reference to start of exon sequence. Case 2.
            exonStart = intronBound.pos + 1
            seqIndex = exonStart - transcript_start + 1
            # refAdd = genome.sequence(
            #     {'chr': chrom, 'start': exonStart - abs(d), 'stop': exonStart - 1}, one_based=True)
            refAdd = genome.get_seq(chrom, exonStart-abs(d), exonStart-1).seq
            exonSeqs[targetExon] = refAdd + exon
            intronBound.pos += d
        if d > 0:  # Need to subtract from start of exon sequence. Case 4
            exonSeqs[targetExon] = exon[d:]
            intronBound.pos += d

        # Modify exon string
        exonCIGARs[targetExon] = editExonCIGAR(exonCIGARs[targetExon], 0, -d)
        intronCIGARs[jn_number] += d

    newSeq = str(''.join(exonSeqs))

    # Paste together the new CIGAR string
    newCIGAR = ""
    for i in range(0, len(intronCIGARs)):
        newCIGAR = newCIGAR + exonCIGARs[i] + str(intronCIGARs[i]) + "N"
    newCIGAR = newCIGAR + exonCIGARs[-1]

    if not check_seq_and_cigar_length(newSeq, newCIGAR):
        raise RuntimeError("CIGAR string and sequence are not the same length")

    return newSeq, newCIGAR


def check_CIGAR_validity(CIGAR):
    """ Returns 'True' if CIGAR string passes tests, and "False" otherwise.
        Checks to make sure that
        1) Introns (N) are only ever followed by M operation
        2) Introns never follow an insertion/deletion
    """

    if "-" in CIGAR:
        return False

    operations, counts = splitCIGAR(CIGAR)
    prev = ""
    for op, ct in zip(operations, counts):
        if prev == "N" and op != "M":
            return False
        if prev == "D" and op == "N":
            return False

        prev = op
    return True


def splitCIGAR(CIGAR):
    """ Takes CIGAR string from SAM and splits it into two lists: one with
        capital letters (match operators), and one with the number of bases """

    matchTypes = re.sub('[0-9]', " ", CIGAR).split()
    matchCounts = re.sub('[A-Z]', " ", CIGAR).split()
    matchCounts = [int(i) for i in matchCounts]

    return matchTypes, matchCounts


def editExonCIGAR(exon, side, nBases):
    """ Given an exon CIGAR string, this function adds or subtracts bases from
        the start (side = 0) or end (side = -1)"""

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
        for i in range(0, len(counts)):
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


def init_log_info(sam_fields):
    """ Initialize a log_info struct to track the error types and correction
        status for a single transcript """

    transcript_id = sam_fields[0]
    flag = int(sam_fields[1])
    chrom = sam_fields[2]

    logInfo = Struct()
    logInfo.TranscriptID = transcript_id

    # Check if read is unmapped, uniquely mapped, or multimapped
    if chrom == "*" or flag == 4:
        logInfo.Mapping = "unmapped"
    elif flag > 16:
        logInfo.Mapping = "non-primary"
    else:
        logInfo.Mapping = "primary"

    logInfo.corrected_deletions = "NA"
    logInfo.uncorrected_deletions = "NA"
    logInfo.variant_deletions = "NA"
    logInfo.corrected_insertions = "NA"
    logInfo.uncorrected_insertions = "NA"
    logInfo.variant_insertions = "NA"
    logInfo.corrected_mismatches = "NA"
    logInfo.uncorrected_mismatches = "NA"
    logInfo.corrected_NC_SJs = "NA"
    logInfo.uncorrected_NC_SJs = "NA"

    return logInfo


def create_log_string(logInfo):
    """ Create string version of logInfo struct """
    log_strings = [str(x) for x in [logInfo.TranscriptID, logInfo.Mapping,
                                    logInfo.corrected_deletions,
                                    logInfo.uncorrected_deletions,
                                    logInfo.variant_deletions,
                                    logInfo.corrected_insertions,
                                    logInfo.uncorrected_insertions,
                                    logInfo.variant_insertions,
                                    logInfo.corrected_mismatches,
                                    logInfo.uncorrected_mismatches,
                                    logInfo.corrected_NC_SJs,
                                    logInfo.uncorrected_NC_SJs]]
    return "\t".join(log_strings)


def write_to_transcript_log(logInfo, tL):
    """ Write a transcript log entry to output """

    log_str = create_log_string(logInfo)
    tL.write(log_str + "\n")
    return


def dryRun(sam, options, outfiles):
    """Records all mismatches, insertions, and deletions in the transcripts,
       but does not correct them """

    print("Dry run mode: Identifying indels and mismatches.........")
    tL = outfiles.log
    eL = outfiles.TElog
    genome = Fasta(options.refGenome)
    sjAnnot = set()

    for line in sam:
        # Set up log entry for transcript
        transcript, logInfo = transcript_init(line, genome, sjAnnot)

        if transcript == None:
            write_to_transcript_log(logInfo, tL)
            continue

        # For primary alignments, modify logInfo
        logInfo.uncorrected_deletions = 0
        logInfo.uncorrected_insertions = 0
        logInfo.uncorrected_mismatches = 0

        # Iterate over CIGAR to catalogue indels and mismatches
        seqPos = 0
        genomePos = transcript.POS
        mergeOperations, mergeCounts = transcript.mergeMDwithCIGAR()

        for op, ct in zip(mergeOperations, mergeCounts):
            if op == "M":
                seqPos += ct
                genomePos += ct

            if op == "X":
                ID = "_".join([transcript.CHROM, str(genomePos)])
                eL.write("\t".join([transcript.QNAME, ID, "Mismatch",
                                    str(ct), "Uncorrected", "DryRun"]) + "\n")
                logInfo.uncorrected_mismatches += 1
                seqPos += ct
                genomePos += ct

            if op == "D":
                ID = "_".join([transcript.CHROM, str(genomePos - 1),
                               str(genomePos + ct - 1)])
                eL.write("\t".join([transcript.QNAME, ID, "Deletion",
                                    str(ct), "Uncorrected", "DryRun"]) + "\n")
                logInfo.uncorrected_deletions += 1
                genomePos += ct

            if op == "I":
                ID = "_".join([transcript.CHROM, str(genomePos - 1),
                               str(genomePos + ct - 1)])
                eL.write("\t".join([transcript.QNAME, ID, "Insertion",
                                    str(ct), "Uncorrected", "DryRun"]) + "\n")
                logInfo.uncorrected_insertions += 1
                seqPos += ct

            if op == "S":
                seqPos += ct

            if op in ["N", "H"]:
                genomePos += ct
        write_to_transcript_log(logInfo, tL)
    return


if __name__ == '__main__':
    #pr = cProfile.Profile()
    # pr.enable()
    main()
    # pr.disable()
    # pr.print_stats(sort='cumtime')
