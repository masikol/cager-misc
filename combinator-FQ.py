#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

__version__ = "1.0.a"
# Year, month, day
__last_update_date__ = "2020-03-02"

# |===== Check python interpreter version =====|

import sys

if sys.version_info.major < 3:
    print( "\nYour python interpreter version is " + "%d.%d" % (sys.version_info.major,
        sys.version_info.minor) )
    print("   Please, use Python 3.\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith("win"):
        raw_input("Press ENTER to exit:")
    # end if
    sys.exit(1)
# end if


def platf_depend_exit(exit_code):
    """
    Function asks to press ENTER press on Windows
        and exits after that.

    :type exit_code: int;
    """
    if sys.platform.startswith("win"):
        input("Press ENTER to exit:")
    # end if
    sys.exit(exit_code)
# end def platf_depend_exit


def err_fmt(text=""):
    """Function for configuring error messages"""
    return "\n  \a!! - ERROR: " + text + '\n'
# end def err_fmt


if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print("combinator-FQ. Version {}. {} edition.".format(__version__, __last_update_date__))
    print("Script screens ends of contigs for identity in order to facilitate further manual scaffolding.")
    print("\nIt generates output table 'combinator_output_FQ.tsv' of following format:")
    print("#  Contig name  Length  Coverage  GC(%)  Multiplicity  Annotation  Start      End")
    print("1  NODE_2...    304356  93.7155   34.12  1.2                       S_NODE_22  S_NODE_26")
    platf_depend_exit(0)
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if


import getopt

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvp:k:i:a:d:",
        ["help", "version", "prefix=", "k-mer-len=", "mink=", "maxk=", "indir="])
except getopt.GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

import re
import os
import glob

prefix = ""
mink = 21
maxk = 127

is_fasta = lambda f: not re.search(r"(m)?f(ast)?a(\.gz)?$", f) is None


if len(args) == 0:

    cont_fpaths = glob.glob( os.path.join(os.getcwd(), "*contigs.f*") )
    cont_fpaths = tuple(filter( is_fasta, cont_fpaths ))

    if len(cont_fpaths) != 0:
        contigs_fpath = cont_fpaths[0]
        print("File '{}' is found and will be processed.".format(contigs_fpath))
        error = True
        while error:
            reply = input("""Press ENTER to process this file
  or enter 'q' to quit:>> """)
            if reply == "":
                error = False
            elif reply == 'q':
                sys.exit(0)
            else:
                print("Invalid reply: '{}'\n".format(reply))
            # end if
        # end while

    else:
        print("No files considered as input fot combinator-FQ is are found.")
        print("Exitting...")
        platf_depend_exit(1)
    # end if

else:
    contigs_fpath = args[0]
    if not os.path.exists(contigs_fpath):
        print("File '{}' does not exist.".format(contigs_fpath))
        platf_depend_exit(1)
    # end if
# end if

for opt, arg in opts:

    if opt in ("-k", "--k-mer-len"):

        try:
            arg = int(arg)
            if arg <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\nLength of k-mer must be positive integer number.")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        else:
            mink, maxk = arg, arg
        # end try

    elif opt in ("-i", "--mink"):

        try:
            mink = int(arg)
            if mink <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\nMinimum length of k-mer must be positive integer number.")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-a", "--maxk"):

        try:
            maxk = int(arg)
            if maxk <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\nMaximum length of k-mer must be positive integer number.")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-p", "--prefix"):
        prefix = arg
    # end if
# end for


if mink > maxk:
    print(err_fmt("Minimum length of k-mer is greater than maximum length of k-mer."))
    print("Values specified by you:")
    print("Minimum length of k-mer: {}.".format(mink))
    print("Maximum length of k-mer: {}.".format(maxk))
    platf_depend_exit(1)
# end if


# Dictionary maps complementary bases according to IUPAC:
compl_dict = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'U': 'A', 'N': 'N',
    'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
    'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
    'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
    'h': 'd', 'v': 'b', 'u': 'a', 'n': 'n'
}

# Function returns complement "comrade" of given base
get_compl_base = lambda base: compl_dict[base]

# Function returns reverse-complement "comrade" of passed DNA sequence
rc = lambda seq: "".join(map( get_compl_base, reversed(seq) ))


from gzip import open as open_as_gzip

is_gzipped = lambda file: True if file.endswith(".gz") else False

# For opening plain text and gzipped files
OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode("utf-8").strip()  # format gzipped line
)


print("\nFile '{}' if processing...".format(contigs_fpath))

open_func = OPEN_FUNCS[ is_gzipped(contigs_fpath) ]
fmt_func = FORMATTING_FUNCS[ is_gzipped(contigs_fpath) ]

contigs_names = list() # list for names of contigs
contigs_seqs = list() # list for sequences of contigs
seq = ""

with open_func(contigs_fpath) as contigs_file:

    for line in contigs_file: # line by line
        line = fmt_func(line)
        if line[0].startswith('>'):
            contigs_names.append(line[1:])
            contigs_seqs.append(seq)
            seq = ""
        else:
            seq += line
        # end if
    # end for
# end with

contigs_seqs.append(seq) # append last sequence to list
contigs_seqs = contigs_seqs[1:] # remove first empty string

# Pattern that matches ID of seqeunce in FASTA file generated by SPAdes
spades_patt = r"^NODE_[0-9]+"
# Pattern that matches ID of seqeunce in FASTA file generated by A5
a5_patt = r"^scaffold_[0-9]+"


voc_ends = dict() # dictionary for collecting contigs data
# Variables for collecting statistics
total_length = 0
agv_coverage = 0
min_coverage = float("inf")
max_coverage = 0

# Readable indices
FULL_NAME, NAME, LEN, COV, GC, START, RC_START, END, RC_END, START_MATCH, END_MATCH = range(11)

# Iterate over contigs and form voc_ends dictionary
for i, contig_name in enumerate(contigs_names):

    # SPAdes and A5 contigs should be processed in different way
    if not re.search(spades_patt, contig_name) is None:
        # Parse fasta header:
        cov = str(contig_name.split('_')[5]) # get coverage
        name = 'NODE_' + contig_name.split('_')[1] # get name in 'NODE_<NUMBER>' format
    elif not re.search(a5_patt, contig_name) is None:
        # Parse fasta header:
        cov = '-' # A5 does not give us coverage
        name = contig_name # use full header as name
    else:
        print("Format of fasta header not recognized:")
        print(" '{}'".format(contig_name))
        print("Please, use SPAdes or A5 assembly.")
        platf_depend_exit(1)
    # end if

    contig_len = len(contigs_seqs[i])

    # Calculating GC-content
    gc_content = 0
    for up_base, low_base in zip(('G', 'C', 'S'),('g', 'c', 's')):
        gc_content += contigs_seqs[i].count(up_base) + contigs_seqs[i].count(low_base)
    # end for

    gcContent = round((gc_content / len(contigs_seqs[i]) * 100), 2)


    # Some common info
    voc_ends.setdefault(i, []).append(contig_name)              # FULL_NAME
    voc_ends.setdefault(i, []).append(name)                         # NAME
    voc_ends.setdefault(i, []).append(contig_len)                   # LEN
    voc_ends.setdefault(i, []).append(cov)                          # COV
    voc_ends.setdefault(i, []).append(gcContent)                    # GC

    # Ends of contigs:
    voc_ends.setdefault(i, []).append(contigs_seqs[i][:maxk])       # START
    voc_ends.setdefault(i, []).append(rc(contigs_seqs[i][:maxk]))   # RC_START
    voc_ends.setdefault(i, []).append(contigs_seqs[i][-maxk:])      # END
    voc_ends.setdefault(i, []).append(rc(contigs_seqs[i][-maxk:]))  # RC_END

    # Ends' matches
    voc_ends.setdefault(i, []).append(list())                       # START_MATCH
    voc_ends.setdefault(i, []).append(list())                       # END_MATCH

    # Collecting statistics
    flt_cov = float(cov)

    total_length += contig_len
    agv_coverage += flt_cov

    min_coverage = min(min_coverage, flt_cov)
    max_coverage = max(max_coverage, flt_cov)
# end for

del contigs_seqs

# Total number of contigs
N = len(contigs_names)
# Expected length of the genome will be less than sum of contig lengths.
# Length of overlapping sticky ends will be subtracted from total_length
exp_genome_len = total_length


def find_overlap(seq1, seq2, mink, maxk):
    """
    Function searches for identity between seq1 and seq2.
    Function regards overlap of length [mink, maxk].

    :param seq1: one sequence;
    :type seq1: str;
    :param seq2: another sequence;
    :type seq2: str;
    :param mink: minimun overlap;
    :type mink: int;
    :param maxk: maximum overlap;
    :type maxk: int;

    Returns 0 if overlap is less than 'mink' and
      length of overlap (which is <= maxk) otherwise.
    """

    i = mink
    overlap = 0

    while seq1[:i] == seq2[:i] and i <= maxk:
        overlap += 1
        i += 1
    # end while

    return 0 if overlap == 0 else (mink + overlap - 1)
# end def find_overlap


print('\n' + '=' * 45)

# Find matching ends of contigs
print("Searching for matching ends of contigs...")
log = "\nMatching ends of contigs:\n"

for i in range(len(contigs_names)):

    if voc_ends[i][LEN] <= mink:
        break
    # end if

    print('\n' + '-'*45 + '\n')
    log += '\n' + '-'*45 + '\n'
    print(voc_ends[i][NAME])
    log += voc_ends[i][NAME] + '\n'

    # === Compare start to end of current contig ===
    # Such contigs probably are adjacent -- write this match to table
    overlap = find_overlap(voc_ends[i][START],
        voc_ends[i][END], mink, maxk)
    if overlap != 0 and voc_ends[i][LEN] != overlap:
        voc_ends[i][START_MATCH].append("[Circle; ovl={}]".format(overlap))
        voc_ends[i][END_MATCH].append("[Circle; ovl={}]".format(overlap))
        exp_genome_len -= overlap
        str_to_print = "{} is circular with overlap of {} b.p.".format(voc_ends[i][FULL_NAME],
            overlap)
        print(str_to_print)
        log += str_to_print + '\n'
    # end if

    # === Compare start to rc-end of current contig ===
    # Contigs probably are not adjacent -- do not write to table
    overlap = find_overlap(voc_ends[i][START],
        voc_ends[i][RC_END], mink, maxk)
    if overlap != 0:
        str_to_print = "Start of {} is identical to it's own rc-end with overlap of {} b.p.".format(voc_ends[i][NAME],
            overlap)
        print(str_to_print)
        log += str_to_print + '\n'
    # end if


    # |=== Compare i-th contig to contigs from i+1 to N ===|
    #   in order not to compare pairs of contigs more than one time
    for j in range(i+1, N):

        # === Compare i-th start to j-th end ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap(voc_ends[i][START],
            voc_ends[j][END], mink, maxk)
        if overlap != 0:
            voc_ends[i][START_MATCH].append("[Start=End({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][END_MATCH].append("[End=Start({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            exp_genome_len -= overlap
            str_to_print = "End of {} matches start of {} with overlap of {} b.p.".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if

        # === Conpare i-th end to j-th start ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap(voc_ends[i][END],
            voc_ends[j][START], mink, maxk)
        if overlap != 0:
            voc_ends[i][END_MATCH].append("[End=Start({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][START_MATCH].append("[Start=End({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            exp_genome_len -= overlap
            str_to_print = "End of {} matches start of {} with overlap of {} b.p.".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if

        # === Compare i-th start to reverse-complement j-th start ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap(voc_ends[i][START],
            voc_ends[j][RC_START], mink, maxk)
        if overlap != 0:
            voc_ends[i][START_MATCH].append("[Start=RC_Start({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][START_MATCH].append("[Start=RC_Start({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            exp_genome_len -= overlap
            str_to_print = "Start of {} matches rc-start of {} with overlap of {} b.p.".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if

        # === Conpare i-th end to reverse-complement j-th end ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap(voc_ends[i][END],
            voc_ends[j][RC_END], mink, maxk)
        if overlap != 0:
            voc_ends[i][END_MATCH].append("[End=RC_End({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][END_MATCH].append("[End=RC_End({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            exp_genome_len -= overlap
            str_to_print = "End of {} matches rc-end of {} with overlap of {} b.p.".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if

        # === Compare i-th start to j-th start ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap(voc_ends[i][START],
            voc_ends[j][START], mink, maxk)
        if overlap != 0:
            str_to_print = "Start of {} matches start of {} with overlap of {} b.p.".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if

        # === Compare i-th end to j-th end ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap(voc_ends[i][END],
            voc_ends[j][END], mink, maxk)
        if overlap != 0:
            str_to_print = "End of {} matches end of {} with overlap of {} b.p.".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if

        # === Compare i-th start to reverse-complement j-th end ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap(voc_ends[i][START],
            voc_ends[j][RC_END], mink, maxk)
        if overlap != 0:
            str_to_print = "End of {} matches rc-start of {} with overlap of {} b.p.".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if

        # === Compare i-th end to reverse-complement j-th start ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap(voc_ends[i][END],
            voc_ends[j][RC_START], mink, maxk)
        if overlap != 0:
            str_to_print = "End of {} matches rc-start of {} with overlap of {} b.p.".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            print(str_to_print)
            log += str_to_print + '\n'
        # end if
    # end for
# end for


# === Form result table in file "combinator_output_FQ.tsv" ===

# Open output file
with open("{}{}combinator_output_FQ.tsv".format(prefix, '' if prefix == '' else '_'), 'w') as outfile:

    # Write head of the table:
    outfile.write("#\tContig name\tLength\tCoverage\tGC(%)\tMultiplicity\tAnnotation\tStart\tEnd\n")

    rev_LQ = 0 # variable for calculating LQ-coefficient
    for i in range(N):

        # Write name of contig:
        outfile.write(str(i + 1) + '\t')
        outfile.write(voc_ends[i][FULL_NAME] + '\t')

        # Write some common data:
        outfile.write(str(voc_ends[i][LEN]) + '\t')
        outfile.write(str(voc_ends[i][COV]) + '\t')
        outfile.write(str(voc_ends[i][GC]) + '\t')

        # Calculate multiplicity of contig:
        muliplty = round(float(voc_ends[i][COV])/float(voc_ends[0][COV]), 1)
        outfile.write(str(muliplty) + '\t\t') # leave empty field for annotation

        # Write information about discovered adjecency
        outfile.write(' '.join(voc_ends[i][START_MATCH]) + '\t')
        outfile.write(' '.join(voc_ends[i][END_MATCH]) + '\n')

        # Calculate number to subtract from 100 to get LQ-coefficient
        # 1 will be subtracted from 100 each time an end remains unpaired
        rev_LQ += int( len(voc_ends[i][START_MATCH]) == 0 )
        rev_LQ += int( len(voc_ends[i][END_MATCH]) == 0 )
    # end for

    # Summary with some statistics:
    outfile.write('\n' + '=' * 45 + '\n')
    outfile.write("Summary:\n")
    outfile.write("Sum of contig lengths: {} b.p.\n".format(total_length))
    outfile.write("Expected length of the genome: {} b.p.\n".format(exp_genome_len))
    outfile.write("Min coverage: {}.\n".format(min_coverage))
    outfile.write("Max coverage: {}.\n".format(max_coverage))
    outfile.write("Average coverage: {}.\n".format(round(agv_coverage / N, 3)))
    outfile.write("{} contigs were processed.\n".format(N))
    LQ = round(100 - (rev_LQ * 100) / (N * 2), 2)
    outfile.write("LQ-coefficient: {}.\n".format((LQ)))
    outfile.write('=' * 45)

    outfile.write('\n' + log) # append log to end of the file
# end with


print("\nTask is completed!")
platf_depend_exit(0)