#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

__version__ = "1.1.b"
# Year, month, day
__last_update_date__ = "2020-03-09"

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


def platf_depend_exit(exit_code=0):
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

def print_help():
    print("\n    |=== combinator-FQ ===|")
    print("Version {}. {} edition.\n".format(__version__, __last_update_date__))
    print("Script identifies adjacent contigs in order to facilitate further scaffolding.")
    print("""Pair of contigs are considered adjacent if terminus (start/end/reverse-complement start/reverse-complement end)
    of one of them overlaps with terminus of another one.""")
    print("Length of overlap is further referred to as 'k'.""")
    print("Format of input: multi-fasta file containing contigs.")
    print("Script supports contigs assembled by SPAdes and A5.")
    print("""\nScript generates 3 output files:
  1) table containing information about contigs adjecency ('<prefix>_combinator_adjacent_contigs.tsv').
    See https://github.com/masikol/cager-misc/wiki/combinator-FQ for details of this output file format;
  2) file, in which all matchings (not only adjacency-associated) are listed
    ('<prefix>_combinator_full_matching_log.txt');
  3) brief summary ('<prefix>_combinator_summary.txt');""")
    print('='*15 + '\n' + "Options:")
    print("  -h (--help): print help message;\n")
    print("  -v (--version): print version;\n")
    print("""  -i (--mink): minimum k (in b.p.) to consider.
    Integer >0; Default if 21 b.p.\n""")
    print("""  -a (--maxk): maximum k (in b.p.) to consider.
    Integer >0; Default is 127 b.p.\n""")
    print("""  -k (--k-mer-len): exact k (in b.p.). If specified, '-i' and '-a' options are ignored.
    Integer >0; Disabled by default;\n""")
    print("""  -p (--prefix): prefix of output files.
  By default they are named according to name of input file,
  e.g. "contigs_combinartor_summary.txt" etc. for 'contigs.fasta';""")
    print('='*15 + '\n' + "Examples:\n")
    print("  ./combinator-FQ.py contigs.fasta -k 127\n")
    print("  ./combinator-FQ.py another_contigs.fa -i 25 -a 300 -p another_contigs")
    print('\n'+'-'*15)
    print("""If input file is omitted in command, combinator will search for
  first (in alphabetical order) fasta file in working directory and
  ask for conformation to process it:\n""")
    print("  ./combinator-FQ.py -i 25 -a 300")
    print("\nIf 'contigs.fasta' is in working directory, combinator will print this:")
    print("~"*10)
    print("File 'contigs.fasta' is found and will be processed.")
    print("""Press ENTER to process this file
  or enter 'h' to see help massage,
  or enter 'q' to quit:>> """)
    print("~"*10)
# end def print_help

# Firstly check for information-providing flags

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print_help()
    platf_depend_exit()
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit()
# end if


import getopt

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvp:k:i:a:",
        ["help", "version", "prefix=", "k-mer-len=", "mink=", "maxk="])
except getopt.GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

import re
import os
import glob

is_fasta = lambda f: not re.search(r"(m)?f(ast)?a(\.gz)?$", f) is None

# Determine fasta file to process:
if len(args) == 0:
    # If no input file is specified, search in working directory:

    cont_fpaths = sorted(glob.glob( os.path.join(os.getcwd(), "*contigs.f*") ))
    cont_fpaths = tuple(filter( is_fasta, cont_fpaths ))

    if len(cont_fpaths) != 0:
        # If there is a fasta file in working directory -- process it
        contigs_fpath = os.path.abspath(cont_fpaths[0])
        print("File '{}' is found and will be processed.".format(contigs_fpath))
        error = True
        while error:
            reply = input("""Press ENTER to process this file
  or enter 'h' to see help massage,
  or enter 'q' to quit:>> """)
            if reply == "":
                error = False
            elif reply == 'h':
                print_help()
                platf_depend_exit()
            elif reply == 'q':
                sys.exit(0)
            else:
                print("Invalid reply: '{}'\n".format(reply))
            # end if
        # end while

    else:
        # If there are no fasta files in working directory, just print help
        print_help()
        platf_depend_exit()
    # end if

else:
    # Check existance of specified file
    contigs_fpath = os.path.abspath(args[0])
    if not os.path.exists(contigs_fpath):
        print("File '{}' does not exist.".format(contigs_fpath))
        platf_depend_exit(1)
    # end if
# end if


prefix = re.search(r"(.+)\.(m)?f(ast)?a(\.gz)?", os.path.basename(contigs_fpath)).group(1)
mink = 21
maxk = 127

# Parse command-line options

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
            mink, maxk = arg, arg # mink = k and maxk = k
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

# Verify mink and maxk:
if mink > maxk:
    if not "-i" in sys.argv[1:] or not "--mink" in sys.argv[1:]:
        mink = maxk
    elif not "-a" in sys.argv[1:] or not "--maxk" in sys.argv[1:]:
        maxk = mink
    else:
        print(err_fmt("Minimum length of k-mer is greater than maximum length of k-mer."))
        print("Values specified by you:")
        print("Minimum length of k-mer: {}.".format(mink))
        print("Maximum length of k-mer: {}.".format(maxk))
        platf_depend_exit(1)
    # end if
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

print("\n  combinator-FQ (Version {}; {} edition)\n".format(__version__, __last_update_date__))

print("File '{}' is processing...".format(contigs_fpath))

open_func = OPEN_FUNCS[ is_gzipped(contigs_fpath) ]
fmt_func = FORMATTING_FUNCS[ is_gzipped(contigs_fpath) ]

contigs_names = list() # list for names of contigs
contigs_seqs = list() # list for sequences of contigs
seq = "" # temporary variable for collecting sequences

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
contigs_seqs = contigs_seqs[1:] # remove first string: it is empty
del seq

# Pattern that matches ID of seqeunce in FASTA file generated by SPAdes
spades_patt = r"^NODE_[0-9]+"
# Pattern that matches ID of seqeunce in FASTA file generated by A5
a5_patt = r"^scaffold_[0-9]+"


voc_ends = dict() # dictionary for collecting contigs data

# Variables for collecting statistics
total_length = 0
avg_coverage = 0
min_coverage = float("inf")
max_coverage = 0

# Readable indices
FULL_NAME, NAME, LEN, COV, GC, START, RC_START, END, RC_END, START_MATCH, END_MATCH = range(11)


# Iterate over contigs and form voc_ends dictionary
for i, contig_name in enumerate(contigs_names):

    contig_len = len(contigs_seqs[i])
    total_length += contig_len

    # SPAdes and A5 contigs should be processed in different way
    if not re.search(spades_patt, contig_name) is None:
        # Parse fasta header:
        cov = str(contig_name.split('_')[5]) # get coverage

        # Collecting coverage statistics
        flt_cov = float(cov)
        avg_coverage += flt_cov

        min_coverage = min(min_coverage, flt_cov)
        max_coverage = max(max_coverage, flt_cov)
        name = 'NODE_' + contig_name.split('_')[1] # get name in 'NODE_<NUMBER>' format
    elif not re.search(a5_patt, contig_name) is None:
        # Parse fasta header:
        cov = '-' # A5 does not write coverage to fasta headers
        name = contig_name # use full header as name
    else:
        print("Format of fasta header not recognized:")
        print(" '{}'".format(contig_name))
        print("Please, use SPAdes or A5 assembly.")
        platf_depend_exit(1)
    # end if

    # Calculate GC-content
    gc_content = 0
    for up_base, low_base in zip(('G', 'C', 'S'),('g', 'c', 's')):
        gc_content += contigs_seqs[i].count(up_base) + contigs_seqs[i].count(low_base)
    # end for
    gcContent = round((gc_content / len(contigs_seqs[i]) * 100), 2)

    # Some common info
    voc_ends.setdefault(i, []).append(contig_name)                  # FULL_NAME
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
# end for

del contigs_seqs

# Total number of contigs
N = len(contigs_names)

# Expected length of the genome will be less than sum of contig lengths.
# Length of overlapping sticky ends will be subtracted from 'total_length'
exp_genome_len = total_length


def find_overlap_s2s(seq1, seq2, mink, maxk):
    """
    Function searches for identity between starts of seq1 and seq2.
    Function regards overlap of length [mink, maxk].

    :param seq1: one sequence;  :type seq1: str;
    :param seq2: another sequence;  :type seq2: str;
    :param mink: minimun overlap;  :type mink: int;
    :param maxk: maximum overlap;  :type maxk: int;

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
# end def find_overlap_s2s


def find_overlap_e2s(seq1, seq2, mink, maxk):
    """
    Function searches for identity between end of seq1 and start of seq2.
    Function regards overlap of length [mink, maxk].

    :param seq1: one sequence;  :type seq1: str;
    :param seq2: another sequence;  :type seq2: str;
    :param mink: minimun overlap;  :type mink: int;
    :param maxk: maximum overlap;  :type maxk: int;

    Returns 0 if overlap is less than 'mink' and
      length of overlap (which is <= maxk) otherwise.
    """

    for i in range(mink, maxk+1):
        if seq1[-i:] == seq2[:i]:
            return i
        # end if
    # end for

    return 0
# end def find_overlap_e2s


def find_overlap_e2e(seq1, seq2, mink, maxk):
    """
    Function searches for identity between ends of seq1 and seq2.
    Function regards overlap of length [mink, maxk].

    :param seq1: one sequence;  :type seq1: str;
    :param seq2: another sequence;  :type seq2: str;
    :param mink: minimun overlap;  :type mink: int;
    :param maxk: maximum overlap;  :type maxk: int;

    Returns 0 if overlap is less than 'mink' and
      length of overlap (which is <= maxk) otherwise.
    """

    i = mink
    overlap = 0

    while seq1[-i:] == seq2[-i:] and i <= maxk:
        overlap += 1
        i += 1
    # end while

    return 0 if overlap == 0 else (mink + overlap - 1)
# end def find_overlap_e2e

# Find matching ends of contigs
full_log = ""
sys.stdout.write("0/{}".format(N))
sys.stdout.flush()

for i in range(len(contigs_names)):

    # Omit contigs shorter that 'mink'
    if voc_ends[i][LEN] <= mink:
        print("\nContig '{}' is shorter than mink({}). Omitting it...".format(voc_ends[i][NAME], mink))
        full_log += "{}: contig is shorter than mink({}). Omitting it...\n".format(voc_ends[i][NAME], mink)
        continue
    # end if

    # === Compare start to end of current contig ===
    # Such contigs probably are adjacent -- write this match to table
    overlap = find_overlap_e2s(voc_ends[i][END],
        voc_ends[i][START], mink, maxk)
    if overlap != 0 and voc_ends[i][LEN] != overlap:
        voc_ends[i][START_MATCH].append("[Circle; ovl={}]".format(overlap))
        voc_ends[i][END_MATCH].append("[Circle; ovl={}]".format(overlap))
        full_log += "{}: contig is circular with overlap of {} b.p.\n".format(voc_ends[i][NAME],
            overlap)
    # end if

    # === Compare start to rc-end of current contig ===
    # Contigs probably are not adjacent -- do not write to table
    overlap = find_overlap_s2s(voc_ends[i][START],
        voc_ends[i][RC_END], mink, maxk)
    if overlap != 0:
        full_log += "{}: start is identical to it's own rc-end with overlap of {} b.p.\n".format(voc_ends[i][NAME],
            overlap)
    # end if


    # |=== Compare i-th contig to contigs from i+1 to N ===|
    #   in order not to compare pairs of contigs more than one time
    for j in range(i+1, N):

        # === Compare i-th start to j-th end ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap_e2s(voc_ends[j][END],
            voc_ends[i][START], mink, maxk)
        if overlap != 0:
            voc_ends[i][START_MATCH].append("[S=E({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][END_MATCH].append("[E=S({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            full_log += "{}: end matches start of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
            full_log += "{}: start matches end of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
        # end if

        # === Conpare i-th end to j-th start ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap_e2s(voc_ends[i][END],
            voc_ends[j][START], mink, maxk)
        if overlap != 0:
            voc_ends[i][END_MATCH].append("[E=S({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][START_MATCH].append("[S=E({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            full_log += "{}: end matches start of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            full_log += "{}: start matches end of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
        # end if

        # === Compare i-th start to reverse-complement j-th start ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap_e2s(voc_ends[j][RC_START],
            voc_ends[i][START], mink, maxk)
        if overlap != 0:
            voc_ends[i][START_MATCH].append("[S=rc_S({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][START_MATCH].append("[S=rc_S({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            full_log += "{}: start matches rc-start of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            full_log += "{}: start matches rc-start of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
        # end if

        # === Conpare i-th end to reverse-complement j-th end ===
        # Such contigs probably are adjacent -- write this match to table
        overlap = find_overlap_e2s(voc_ends[i][END],
            voc_ends[j][RC_END], mink, maxk)
        if overlap != 0:
            voc_ends[i][END_MATCH].append("[E=rc_E({}); ovl={}]".format(voc_ends[j][NAME],
                overlap))
            voc_ends[j][END_MATCH].append("[E=rc_E({}); ovl={}]".format(voc_ends[i][NAME],
                overlap))
            full_log += "{}: end matches rc-end of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
            full_log += "{}: end matches rc-end of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
        # end if

        # === Compare i-th start to j-th start ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap_s2s(voc_ends[i][START],
            voc_ends[j][START], mink, maxk)
        if overlap != 0:
            full_log += "{}: start matches start of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            full_log += "{}: start matches start of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
        # end if

        # === Compare i-th end to j-th end ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap_e2e(voc_ends[i][END],
            voc_ends[j][END], mink, maxk)
        if overlap != 0:
            full_log += "{}: end matches end of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            full_log += "{}: end matches end of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
        # end if

        # === Compare i-th start to reverse-complement j-th end ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap_s2s(voc_ends[i][START],
            voc_ends[j][RC_END], mink, maxk)
        if overlap != 0:
            full_log += "{}: end matches rc-start of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
            full_log += "{}: start matches rc-end of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
        # end if

        # === Compare i-th end to reverse-complement j-th start ===
        # Contigs probably are not adjacent -- do not write to table
        overlap = find_overlap_e2e(voc_ends[i][END],
            voc_ends[j][RC_START], mink, maxk)
        if overlap != 0:
            full_log += "{}: end matches rc-start of {} with overlap of {} b.p.\n".format(voc_ends[i][NAME],
                voc_ends[j][NAME], overlap)
            full_log += "{}: start matches rc-end of {} with overlap of {} b.p.\n".format(voc_ends[j][NAME],
                voc_ends[i][NAME], overlap)
        # end if
    # end for

    sys.stdout.write("\r{}/{}".format(i+1, N))
    sys.stdout.flush()
# end for


# === Form result table in file "combinator_output_FQ.tsv" ===

# Open output file
with open("{}{}combinator_adjacent_contigs.tsv".format(prefix, '' if prefix == '' else '_'), 'w') as outfile:

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
        if voc_ends[0][COV] != '-':
            muliplty = round(float(voc_ends[i][COV])/float(voc_ends[0][COV]), 1)
        else:
            muliplty = '-'
        # end if

        outfile.write(str(muliplty) + '\t\t') # leave empty field for annotation

        # Write information about discovered adjecency
        for idx, eol in zip((START_MATCH, END_MATCH), ('\t', '\n')):
            str_to_write = ' '.join(voc_ends[i][idx])
            ovls = re.findall(r"; ovl=([0-9]+)", str_to_write)
            if len(ovls) != 0:
                # Subtract maximum overlap in order to find expected genome size
                exp_genome_len -= max(map(lambda x: int(x), ovls))
            else:
                # Calculate number to subtract from 100 to get LQ-coefficient
                # 1 will be subtracted from 100 each time an end remains unpaired
                rev_LQ += 1
            # end if
            outfile.write(str_to_write + eol)
        # end for
    # end for
# end with

# Sort log by name of contig
full_log = '\n'.join(sorted(full_log.splitlines(), key = lambda x: int(x.partition(":")[0].partition('_')[2])))

with open("{}{}combinator_full_matching_log.txt".format(prefix, '' if prefix == '' else '_'), 'w') as outfile:
    outfile.write(full_log) # write full log separate file
# end with

LQ = round(100 - (rev_LQ * 100) / (N * 2), 2)

print('\n')
with open("{}{}combinator_summary_FQ.txt".format(prefix, '' if prefix == '' else '_'), 'w') as outfile:

    outfile.write("File: '{}'\n\n".format(contigs_fpath))

    # Summary with some statistics:
    for print_func in (sys.stdout.write, outfile.write):
        print_func("Summary:\n")
        print_func("{} contigs were processed.\n".format(N))
        print_func("Sum of contig lengths: {} b.p.\n".format(total_length))
        print_func("Expected length of the genome: {} b.p.\n".format(exp_genome_len))
        min_coverage = '-' if min_coverage == float("inf") else min_coverage
        print_func("Min coverage: {}\n".format(min_coverage))
        if max_coverage > 1e-6:
            print_func("Max coverage: {}\n".format(max_coverage))
        else:
            print_func("Max coverage: -\n")
        if avg_coverage > 1e-6:
            print_func("Average coverage: {}\n".format(round(avg_coverage / N, 3)))
        else:
            print_func("Average coverage: -\n")
        # end if
        print_func("LQ-coefficient: {}\n".format((LQ)))
    # end for
# end with
sys.stdout.flush()

platf_depend_exit()