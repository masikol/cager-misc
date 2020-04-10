#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

#---------------------------------------------------------------------------
# Script calculates GC-content of each sequence in *.fasta(.gz) or *.fa(.gz) file(s).
#
# In the end, script prints summary:
# 1. Total length of sequences processed.
# 2. Min, max and mean coverage (if SPAdes assembly file is processed).
#---------------------------------------------------------------------------

__version__ = "1.0.a"
# Year, month, day
__last_update_date__ = "2020-04-10"

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

import os

def print_help():
    print("\n    |=== fasta_GC-content ===|")
    print("Version {}. {} edition.\n".format(__version__, __last_update_date__))
    print("""Script calculates GC-content of each sequence in
  *.fasta or *.fa file(s) and saves this data in file(s) '<NAME>_GC_result.txt'.""")
    print("Also, script can calculate max, min and mean coverage when processing SPAdes's assembly.")

    print("\nUsage:")
    print("  python3 fasta_GC-content.py some_file.fasta another_file.fa")
    print("Following command will process all *.fasta(.gz) and *.fa(.gz) files in the working directory:")
    print("  python3 fasta_GC-content.py")
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

fa_fpaths = list()
import re

is_fasta = lambda f: not re.search(r"f(ast)?a(\.gz)?$", f) is None
valid_options = ("-h", "--help", "-v", "--version", "-q", "--quiet")

for arg in sys.argv[1:]:

    if not arg in valid_options:

        if not os.path.exists(arg):
            print("File '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        # end if

        if not is_fasta(arg):
            print("File '{}' does not like an fasta file".format(arg))
            print("Scrtip understands only '*.fasta(.gz)'' and '*.fa(.gz)' extentions")
            platf_depend_exit(1)
        # end if

        fa_fpaths.append(arg)
# end for

del valid_options


if len(fa_fpaths) == 0:

    # Check *.fasta или *.fa in the working dir
    fa_fpaths = os.listdir(".")
    fa_fpaths = tuple(filter(is_fasta, fa_fpaths))

    # Check if there are any appropriate files:
    if len(fa_fpaths) == 0:
        print("\nThere are no '*.fasta(.gz)' or '*.fa(.gz)' files in the working directory!")
        print_help()
        platf_depend_exit(1)
    # end if
# end if

if "-q" in sys.argv[1:] or "--quiet" in sys.argv[1:]:

    def printq(text=""):
        """Quiet printing"""
        pass
    # end def printq

else:

    def printq(text=""):
        """Normal printing"""
        sys.stdout.write(text + '\n')
    # end def printq

# end if

printq("\nfasta_GC-content.py; Version {}; {} edition;\n".format(__version__, __last_update_date__))
printq('Following files will be processed:')
for i, fpath in enumerate(fa_fpaths):
    printq("{}. '{}'".format(str(i+1), fpath))
printq('-' * 58 + '\n')


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

    printq("========== File: '{}' ===========".format(fpath))
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

    printq("\n{} sequences in file '{}'\n".format(len(heads), fpath))

    # Configure name of outfile with appropriate prefix
    prefix = re.search(r"(.+)\.fa(sta)?(\.gz)?", fpath).group(1)
    outfpath = prefix + ".GC_result.txt"


    with open(outfpath, 'w') as outfile:

        # Calculate GC-content
        for (i, seq), head in zip(enumerate(seqs), heads):

            # Write heads to outfile and console
            outfile.write(head.replace('>', '') + '\t' + '\n')
            printq(head.replace('>', ''))

            g_count = seq.upper().count('G')
            outfile.write("G count:  \t{}\n".format(g_count))
            printq("G count:  \t{}".format(g_count))

            c_count = seq.upper().count('C')
            outfile.write("C count:  \t{}\n".format(c_count))
            printq("C count:  \t{}".format(c_count))

            # According to IUPAC, S is G or C
            s_count = seq.upper().count("S")
            # We wont's disturb a user with this 'S' if s_count == 0
            if s_count != 0:
                outfile.write("S count:  \t{}\n".format(s_count))
                printq("S count:  \t{}".format(s_count))
            # end if

            gc_content = round( (g_count + c_count + s_count) / len(seq) * 100, 2 )
            outfile.write("Length:  \t{} b.p.\n".format(len(seq)))
            printq("Length:  \t{} b.p.".format(len(seq)))
            outfile.write("GC-content:  \t{}%\n".format(gc_content))
            printq("GC-content:  \t{}%".format(gc_content))
            # End of GC-content calculating

            # Calculations for summary
            totalLength += len(seq)

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
        # end for

        # Duplicate printing to console and writing to file
        for print_func in (printq, outfile.write):

            print_func("\n================= Summary ==================\n")
            print_func("Total Length:    \t{:,}\n".format(totalLength).replace(',', ' '))

            if spades:
                print_func("Min coverage:   \t{}\n".format(min_cov))
                print_func("Max coverage:   \t{}\n".format(max_cov))
                print_func("Mean coverage:  \t{}\n".format(round(mean_cov / len(seqs), 3)))
            # end if

            print_func('Sequences processed: \t{}\n'.format('{:,}'.format(len(seqs)).replace(',', ' ')))
            print_func('Lines processed: \t{}\n'.format('{:,}'.format(line_counter).replace(',', ' ')))
            print_func('============================================\n\n')
        # end for
    # end with
# end for

printq("Completed!")
platf_depend_exit()