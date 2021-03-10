#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#---------------------------------------------------------------------------
# Script calculates GC-content of each sequence in *.fasta(.gz) or *.fa(.gz) file(s).
# Script writes it's output to tab-separated file.
#
# For each file, script additionaly prints summary:
# 1. Total length of sequences processed.
# 2. Min, max and mean coverage (if SPAdes assembly file is processed).
#---------------------------------------------------------------------------

__version__ = "1.1.b"
# Year, month, day
__last_update_date__ = "2021-03-10"

# Check python interpreter version

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


def print_help():
    print("\n    |=== fasta-GC-content ===|")
    print("Version {}. {} edition.\n".format(__version__, __last_update_date__))
    print("""Script calculates GC-content of each sequence in
  *.fasta or *.fa file(s) and saves this data in file(s) '<NAME>_GC_result.txt'.""")
    print("Also, script can calculate max, min and mean coverage when processing SPAdes's assembly.")

    print("\nUsage:")
    print("  python3 fasta-GC-content.py some_file.fasta another_file.fa")
    print("Following command will process all *.fasta(.gz) and *.fa(.gz) files in the working directory:")
    print("  python3 fasta-GC-content.py")
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

import re
import os

fa_fpaths = list()

is_fasta = lambda f: not re.search(r"f(ast)?a(\.gz)?$", f) is None
valid_options = ("-h", "--help", "-v", "--version")

for arg in sys.argv[1:]:

    if not arg in valid_options:

        if not os.path.exists(arg):
            print("File '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        # end if

        if not is_fasta(arg):
            print("File '{}' does not look like a fasta file".format(arg))
            print("Script understands only '*.fasta(.gz)'' and '*.fa(.gz)' extentions")
            platf_depend_exit(1)
        # end if

        fa_fpaths.append(arg)
# end for

del valid_options


if len(fa_fpaths) == 0:

    # Check *.fasta *.fa in the working dir
    fa_fpaths = os.listdir(".")
    fa_fpaths = tuple(filter(is_fasta, fa_fpaths))

    # Check if there are any appropriate files:
    if len(fa_fpaths) == 0:
        print_help()
        platf_depend_exit(1)
    # end if
# end if


print("\nfasta-GC-content.py; Version {}; {} edition;\n".format(__version__, __last_update_date__))
print('Following files will be processed:')
for i, fpath in enumerate(fa_fpaths):
    print("{}. '{}'".format(str(i+1), fpath))
print('-' * 25 + '\n')


# Stuff for processing gzipped files

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


# ------------- Start processing ------------

for fpath in fa_fpaths:

    how_to_open = OPEN_FUNCS[ is_gzipped(fpath) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fpath) ]

    print("=== File: '{}' ===".format(fpath))
    seq = ''            # empty string
    heads = []          # list for headings
    seqs = []           # list of sequences
    line_counter = 0    # counter for processed lines

    # Variables for summary
    totalLength = 0
    mean_cov = 0
    min_cov = float('inf')
    max_cov = 0

    with how_to_open(fpath) as infile:

        for line_str in infile:

            line_str = fmt_func(line_str)

            if line_str[0] == '>':             # if line_str is seq's heading
                heads.append(line_str.strip()) # append head to 'heads'
                seqs.append(seq)               # append "accumulated" sequence to 'seqs'
                seq = ''                       # reset seq
            else:
                seq += line_str.strip()        # "accumulate" sequence
            line_counter += 1

        seqs.append(seq)                       # append the last sequence to 'seqs'
    # end with

    # Remove odd empty string (it is the first one)
    seqs.remove('')

    # Configure name of outfile with appropriate prefix
    prefix = re.search(r"(.+)\.fa(sta)?(\.gz)?", fpath).group(1)
    outfpath = prefix + ".GC_result.txt"


    with open(outfpath, 'w') as outfile:

        outfile.write('\t'.join( ("Sequence name",
            "G",
            "C",
            "S (degenerate)",
            "GC (%)",
            "Length (b.p.)",
            "Coverage") ) +
        '\n')

        # Calculate GC-content
        for (i, seq), head in zip(enumerate(seqs), heads):

            # Check if we process SPAdes's assembly
            spades = True
            head_split = head.split('_')
            spades = spades and head_split[0] == '>NODE'
            spades = spades and head_split[2] == 'length'
            spades = spades and head_split[4] == 'cov'

            if spades:
                # SPAdes gracefully bestows coverage upon us
                cov = head_split[5]
                mean_cov += float(cov)

                min_cov = min(min_cov, float(cov))
                max_cov = max(max_cov, float(cov))
            # end if

            g_count = seq.upper().count('G')
            c_count = seq.upper().count('C')
            # According to IUPAC, S is G or C
            s_count = seq.upper().count("S")

            gc_content = round( (g_count + c_count + s_count) / len(seq) * 100, 2 )

            outfile.write('\t'.join( (head[1:],
                str(g_count),
                str(c_count),
                str(s_count),
                str(gc_content),
                str(len(seq)),
                str(cov))) +
            '\n')

            # Calculations for summary
            totalLength += len(seq)
        # end for

        # Separate summary from table body
        outfile.write('\n')

        # Duplicate printing to console and writing to file
        for print_func in (sys.stdout.write, outfile.write):

            print_func("Total length:    \t{:,}\n".format(totalLength).replace(',', ' '))

            if spades:
                print_func("Min coverage:   \t{}\n".format(min_cov))
                print_func("Max coverage:   \t{}\n".format(max_cov))
                print_func("Mean coverage:  \t{}\n".format(round(mean_cov / len(seqs), 3)))
            # end if

            print_func('Sequences processed: \t{}\n'.format('{:,}'.format(len(seqs)).replace(',', ' ')))
            print_func('Lines processed: \t{}\n'.format('{:,}'.format(line_counter).replace(',', ' ')))
        # end for
        print('=' * 30 + '\n')
    # end with
# end for

print("Completed!")
platf_depend_exit()