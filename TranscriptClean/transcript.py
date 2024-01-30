# This file contains the Transcript2 class and associated methods for the
# TranscriptClean program

from spliceJunction import SpliceJunction
import re
import itertools


class Transcript:
    def __init__(self, samFields, genome, spliceAnnot):

        # These eleven attributes are initialized directly from the input
        # SAM entry and are mandatory
        self.QNAME = samFields[0]
        self.FLAG = int(samFields[1])
        self.CHROM = samFields[2]
        self.POS = int(samFields[3])
        self.MAPQ = samFields[4]
        self.CIGAR = str(samFields[5])
        self.RNEXT = samFields[6]
        self.PNEXT = samFields[7]
        self.TLEN = samFields[8]
        self.SEQ = samFields[9]
        self.QUAL = "*"

        # If the sam entry contains additional optional fields, process them
        self.NM = ""
        self.MD = ""
        self.jM = ""
        self.jI = ""
        otherFields = []

        for field in samFields[11:len(samFields)]:
            if field.startswith("NM"):
                self.NM = field
            elif field.startswith("MD"):
                self.MD = field
            elif field.startswith("jM"):
                self.jM = field
            elif field.startswith("jI"):
                self.jI = field
            else:
                otherFields.append(field)

        # If the NM and MD tags are None, it means there was a reference genome
        # problem somewhere in the read. Consider such reads unmapped.
        if self.NM == "" or self.MD == "":
            self.NM, self.MD = self.getNMandMDFlags(genome)

        # These attributes are set by parsing the inputs
        self.strand = "+"
        if int(self.FLAG) == 16 or int(self.FLAG) == 272:
            self.strand = "-"

        # Get intron locations from the CIGAR string if not included already
        if (self.jI == ""):
            self.jI = self.compute_jI()

        if "N" in self.CIGAR:
            # Create an object for each splice junction
            self.spliceJunctions, self.isCanonical, self.allJnsAnnotated = \
                self.parseSpliceJunctions(genome, spliceAnnot)
        else:
            self.spliceJunctions = []
            self.isCanonical = True
            self.allJnsAnnotated = True

        # Get annotation status of each junction
        if (self.jM == ""):
            self.jM, self.jI = self.get_jM_jI_tags_from_sjs()

        self.otherFields = "\t".join(otherFields)

    def recheckJnsAnnotated(self):
        """ Check the splice motif of each splice junction to determine
            whether the transcript overall is annotated """
        for jn in self.spliceJunctions:
            if int(jn.motif_code) < 20:
                self.allJnsAnnotated = False
                return False
        self.allJnsAnnotated = True
        return True

    def recheckCanonical(self):
        """ Check each splice junction. If one or more junctions are
            noncanonical, then so is the transcript. """
        for jn in self.spliceJunctions:
            if jn.isCanonical == False:
                self.isCanonical = False
                return False
        self.isCanonical = True
        return True

    def compute_transcript_end(self):
        """ Given the start position and CIGAR string of a mapped SAM transcript,
            compute the end position in the reference genome.

            Args:
                start: The start position of the transcript with respect to the
                forward strand

                cigar: SAM CIGAR string describing match operations to the reference
                genome

            Returns:
                end position of the transcript.
        """
        end = self.POS

        ops, counts = self.splitCIGAR()
        for op, ct in zip(ops, counts):
            if op in ["M", "N", "D"]:
                end += ct

        return end - 1

    def splitCIGAR(self):
        """ Takes CIGAR string from SAM and splits it into two lists:
            one with capital letters (match operators), and one with
            the number of bases that each operation applies to. """

        #alignTypes = re.sub('[0-9]', " ", self.CIGAR).split()
        #counts = re.sub('[A-Z]', " ", self.CIGAR).split()
        #counts = [int(i) for i in counts]

        return splitCIGARstr(self.CIGAR)  # alignTypes, counts

    def splitMD(self):
        """ Takes MD tag and splits into two lists:
            one with capital letters (match operators), and one with
            the number of bases that each operation applies to. """

        MD = str(self.MD).split(":")[2]
        operations = []

        # Split MD string where type changes.
        # Digits are separated from base changes.
        # Deletions (with ^) are captured together.
        counts = ["".join(x)
                  for _, x in itertools.groupby(MD, key=str.isdigit)]

        # Get operations
        for i in range(0, len(counts)):
            curr = counts[i]
            try:
                counts[i] = int(curr)
                operations.append("M")
            except ValueError:
                # Handle deletion
                if curr.startswith("^"):
                    operations.append("D")
                    counts[i] = len(counts[i]) - 1
                else:
                    operations.append("X")
                    counts[i] = len(counts[i])

        return operations, counts

    def mergeMDwithCIGAR(self):
        """ Takes the MD and CIGAR strings, and combines them into a unified
            structure that encodes all possible operations w.r.t the reference:
            match, mismatch, deletion, insertion, hard clipping,
            and soft clipping. """

        mergeCounts = []
        mergeOperations = []

        cigarOperation, cigarCount = self.splitCIGAR()
        mdOperation, mdCount = self.splitMD()

        mdIndex = 0
        cigarIndex = 0

        while mdIndex < len(mdOperation) or cigarIndex < len(cigarOperation):

            # If the current CIGAR operation is S, H, N, or I, add that to the
            # output. The MD tag doesn't have these
            if cigarOperation[cigarIndex] in ("H", "S", "I", "N"):
                mergeOperations.append(cigarOperation[cigarIndex])
                mergeCounts.append(cigarCount[cigarIndex])
                cigarIndex += 1

            # Otherwise, we need to compare the current CIGAR and MD operations.
            # Select the "shorter" operation and add it to the results.
            # Subtract away the same number of bases from the competing entry.
            else:
                if cigarCount[cigarIndex] < mdCount[mdIndex]:
                    # If the CIGAR string lists fewer matched bases than MD,
                    # it means the CIGAR has had an insertion not listed in MD
                    mdCount[mdIndex] = mdCount[mdIndex] - \
                        cigarCount[cigarIndex]
                    mergeOperations.append(cigarOperation[cigarIndex])
                    mergeCounts.append(cigarCount[cigarIndex])
                    cigarIndex += 1

                elif cigarCount[cigarIndex] > mdCount[mdIndex]:
                    # If the CIGAR string lists more matched bases than MD,
                    # it means that MD has a mismatch not listed in CIGAR
                    cigarCount[cigarIndex] = cigarCount[cigarIndex] - \
                        mdCount[mdIndex]
                    mergeOperations.append(mdOperation[mdIndex])
                    mergeCounts.append(mdCount[mdIndex])
                    mdIndex += 1

                # For cases where both MD and CIGAR specify the same match type,
                # add to the result and advance to next position in lists
                else:
                    mergeOperations.append(mdOperation[mdIndex])
                    mergeCounts.append(mdCount[mdIndex])
                    mdIndex += 1
                    cigarIndex += 1

        return mergeOperations, mergeCounts

    def parseSpliceJunctions(self, genome, spliceAnnot):
        """ Takes the splice junction information from the SAM input and
            creates a SpliceJunction object for each junction."""

        intronBounds = ((self.jI).split(":")[-1]).split(",")[1:]

        count = 0
        jnNum = 0
        jnObjects = []
        canonical = True
        annotated = True
        while count < len(intronBounds):
            start = int(intronBounds[count])
            end = int(intronBounds[count + 1])
            sj = SpliceJunction(self.QNAME, jnNum, self.CHROM, start, end,
                                self.strand, genome, spliceAnnot)
            jnObjects.append(sj)

            # Check if junction is canonical or not, as well as whether it is
            # annotated.
            if sj.isCanonical == False:  # self.isCanonical = False
                canonical = False
            if int(sj.motif_code) < 20:
                annotated = False
            count += 2
            jnNum += 1

        return jnObjects, canonical, annotated

    def printableSAM(self):
        """ Returns a SAM-formatted string representation of the transcript"""
        fields = [self.QNAME, self.FLAG, self.CHROM, self.POS, self.MAPQ, self.CIGAR,
                  self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, self.QUAL, self.otherFields,
                  self.NM, self.MD, self.jM, self.jI]

        final_fields = []
        for field in fields:
            if field != "" and field != None:
                final_fields.append(field)

        return "\t".join([str(x) for x in final_fields]).strip()

    def printableFa(self):
        """ Returns a fasta-formatted string representation of the transcript """
        fastaID = ">" + self.QNAME
        strand = self.strand
        seq = str(self.SEQ)

        if strand == "-":  # Need to reverse-complement the sequence
            seq = reverseComplement(seq)

        # Split seq into 80-character segments
        fastaSeq = [seq[i:i+80] for i in range(0, len(seq), 80)]
        return fastaID + "\n" + "\n".join(fastaSeq)

    def getAllIntronBounds(self):
        """ Return all intron bound objects belonging to this transcript """

        result = []
        for jn in self.spliceJunctions:
            b = jn.bounds
            result.append(b[0])
            result.append(b[1])
        return result

    def getNMandMDFlags(self, genome):
        """ This function uses the transcript sequence, its CIGAR string,
            and the reference genome to create NM and MD sam flags."""
        NM = 0
        MD = "MD:Z:"
        MVal = 0
        seqPos = 0
        genomePos = self.POS

        operations, counts = self.splitCIGAR()

        tot = 0
        for op, ct in zip(operations, counts):
            if op in ["M", "I", "S"]:
                tot += ct

        for op, ct in zip(operations, counts):
            if op == "M":
                for i in range(0, ct):
                    currBase = self.SEQ[seqPos]
                    refBase = genome.get_seq(self.CHROM, genomePos, genomePos).seq
                    # refBase = genome.sequence({'chr': self.CHROM, 'start': genomePos,
                    #                            'stop': genomePos}, one_based=True)
                    if refBase == "":
                        return None, None

                    # In the event of a mismatch
                    if currBase.upper() != refBase.upper():
                        # End any match we have going and add the mismatch
                        MD = MD + str(MVal)
                        MVal = 0
                        MD = MD + str(refBase)
                        NM += 1
                    # Bases match
                    else:
                        MVal += 1
                    # Either way, advance forwards in the sequence and genome
                    seqPos += 1
                    genomePos += 1
            if op == "D":
                # End any match we have going and add the missing reference bases
                MD = MD + str(MVal)
                MVal = 0
                refBases = genome.get_seq(self.CHROM, genomePos,
                                          genomePos+ct-1).seq
                # refBases = genome.sequence({'chr': self.CHROM, 'start': genomePos,
                #                             'stop': genomePos + ct - 1}, one_based=True)
                if refBases == "":
                    return None, None
                MD = MD + "^" + str(refBases)
                NM += ct
                genomePos += ct
            # For insertions and soft clips, we move on without adding to the MD
            if op in ["I", "S"]:
                seqPos += ct
                if op == "I":
                    NM += ct
            if op in ["N", "H"]:
                genomePos += ct

        if MVal > 0:
            MD = MD + str(MVal)
        return "NM:i:" + str(NM), MD

    def get_jM_jI_tags_from_sjs(self):
        """ Create jM and jI tags by traversing the splice junction strings """

        jM = ["jM:B:c"]
        jI = ["jI:B:i"]

        for sj in self.spliceJunctions:
            intron_start = sj.bounds[0].pos
            intron_end = sj.bounds[1].pos
            motif_code = sj.motif_code
            jM.append(str(motif_code))
            jI.append(str(intron_start))
            jI.append(str(intron_end))

        # If the transcript has no introns, we need to add -1 to the tags
        if len(jM) == len(jI) == 1:
            jM.append("-1")
            jI.append("-1")

        jMstr = ",".join(jM)
        jIstr = ",".join(jI)

        return jMstr, jIstr

    def compute_jI(self):
        """ Use the CIGAR string to compute where the introns are """

        operations, counts = self.splitCIGAR()
        jI = ["jI:B:i"]
        genomePos = self.POS

        # Iterate over operations
        for op, ct in zip(operations, counts):
            if op == "N":
                # This is an intron
                intronStart = genomePos
                intronEnd = genomePos + ct - 1
                jI.append(str(intronStart))
                jI.append(str(intronEnd))

            if op not in ["S", "I"]:
                genomePos += ct

        # If the transcript has no introns, we need to add -1 to the tags
        if len(jI) == 1:
            jI.append("-1")

        jIstr = ",".join(jI)

        return jIstr

    def getjMandjITags(self, genome, spliceAnnot):
        """ If the input sam file doesn't have the custom STARlong-derived jM
            and jI tags, we need to compute them. This is done by stepping
            through the CIGAR string and sequence. When an intron (N) is
            encountered, we check the first two bases and last two bases of
            the intron in the genome sequence to detemine whether they are
            canonical. We also record the start and end position of the intron. """

        seq = self.SEQ
        operations, counts = self.splitCIGAR()

        jM = ["jM:B:c"]
        jI = ["jI:B:i"]

        genomePos = self.POS

        # Iterate over operations
        for op, ct in zip(operations, counts):
            if op == "N":
                # This is an intron
                intronStart = genomePos
                # startBases = genome.sequence({'chr': self.CHROM,
                #                               'start': genomePos,
                #                               'stop': genomePos + 1},
                #                              one_based=True)
                startBases = genome.get_seq(self.CHROM, genomePos,
                                            genomePos+1).seq
                intronEnd = genomePos + ct - 1
                # endBases = genome.sequence({'chr': self.CHROM,
                #                             'start': intronEnd - 1,
                #                             'stop': intronEnd}, one_based=True)
                endBases = genome.get_seq(self.CHROM, intronEnd-1, intron_end).seq

                # Check if junction is annotated
                if self.strand == "+":
                    type1 = "donor"
                    type2 = "acceptor"
                elif self.strand == "-":
                    type1 = "acceptor"
                    type2 = "donor"

                if ("_".join([self.CHROM, str(intronStart), self.strand, type1])) in spliceAnnot and \
                   ("_".join([self.CHROM, str(intronEnd), self.strand, type2])) in spliceAnnot:
                    motifCode = 20 + getSJMotifCode(startBases, endBases)
                else:
                    motifCode = getSJMotifCode(startBases, endBases)

                jM.append(str(motifCode))
                jI.append(str(intronStart))
                jI.append(str(intronEnd))

            if op not in ["S", "I"]:
                genomePos += ct

        # If the transcript has no introns, we need to add -1 to the tags
        if len(jM) == len(jI) == 1:
            jM.append("-1")
            jI.append("-1")

        jMstr = ",".join(jM)
        jIstr = ",".join(jI)

        return jMstr, jIstr

    def base_wise_CIGAR(self):
        """ Create an extended version of the CIGAR string with an operation
            per base.
            Example: 3M10N2D5M becomes
                     MMMNNNNNNNNNNDDMMMMM

        """
        base_wise_str = ""
        ops, counts = self.splitCIGAR()
        for op, ct in zip(ops, counts):
            base_wise_str += op*ct
        return base_wise_str

    def fetch_region_sequence(self, chromosome, start, end):
        """ Walks the SAM sequence to return the bases in the specified region.
            Returns None if the sequence is not available, i.e. because the
            region does not overlap with the transcript, or because it is in
            an intron. Supplied coordinates should be 1-based. """

        # Check whether the supplied region is located even remotely near the
        # transcript
        if chromosome != self.CHROM:
            return None
        if not(start >= self.POS and end <= self.compute_transcript_end()):
            return None

        # Walk transcript sequence using CIGAR string
        positions = range(start, end + 1)
        seq = self.SEQ
        seq_pos = 0
        bases = ""

        genome_pos = self.POS
        while genome_pos <= end:
            for op in self.base_wise_CIGAR():
                if op == "M":
                    if (genome_pos in positions):
                        bases += seq[seq_pos]
                    genome_pos += 1
                    seq_pos += 1
                # Advance in genome sequence but not in transcript sequence
                if op in ["D", "N", "H"]:
                    if (genome_pos in positions) and op == "D":
                        bases += "-"
                    genome_pos += 1
                # Advance in transcript sequence but not genome sequence
                if op in ["S", "I"]:
                    if (genome_pos in positions) and op == "I":
                        bases += seq[seq_pos]
                    seq_pos += 1

        if bases == "":
            return None

        # If the transcript is on the reverse strand, reverse-complement the
        # sequence before returning it
        if self.strand == "-":
            bases = reverseComplement(bases)
        return bases


def getSJMotifCode(startBases, endBases):
    """ Determines which STAR-style splice junction code applies to a splice motif """

    motif = (startBases + endBases).upper()

    if motif == "GTAG":
        return 1
    elif motif == "CTAC":
        return 2
    elif motif == "GCAG":
        return 3
    elif motif == "CTGC":
        return 4
    elif motif == "ATAC":
        return 5
    elif motif == "GTAT":
        return 6
    else:
        return 0


def reverseComplement(seq):
    """ Returns the reverse complement of a DNA sequence,
        retaining the case of each letter"""
    complement = ""

    for base in seq:
        if base == "A":
            complement += "T"
        elif base == "T":
            complement += "A"
        elif base == "G":
            complement += "C"
        elif base == "C":
            complement += "G"
        elif base == "N":
            complement += "N"
        elif base == "a":
            complement += "t"
        elif base == "t":
            complement += "a"
        elif base == "g":
            complement += "c"
        elif base == "c":
            complement += "g"
        elif base == "n":
            complement += "n"
        elif base == "*":
            complement += "*"
        else:
            complement += base
            print(
                "Warning: reverse complement function encountered unknown base " + "'" + base + "'")

    reverseComplement = complement[::-1]

    return reverseComplement


def splitCIGARstr(CIGAR):
    """ Takes CIGAR string from SAM and splits it into two lists:
        one with capital letters (match operators), and one with
        the number of bases that each operation applies to. """

    alignTypes = re.sub('[0-9]', " ", CIGAR).split()
    counts = re.sub('[A-Z]', " ", CIGAR).split()
    counts = [int(i) for i in counts]

    return alignTypes, counts


def check_seq_and_cigar_length(seq, cigar):
    """This function computes the sequence and CIGAR length, then compares
       them to see if they are the same. Returns True if yes, False if not.
    """
    seq_len = len(seq)
    ops, counts = splitCIGARstr(cigar)
    cigar_len = 0
    for op, ct in zip(ops, counts):
        if op in ["M", "I", "S"]:
            cigar_len += int(ct)

    if seq_len == cigar_len:
        return True
    else:
        return False
