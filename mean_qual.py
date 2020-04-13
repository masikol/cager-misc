#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------
# This script calculates total mean quality of reads in fasrq file(s).
# Script writes it's output to tab-separated file.
# ----------------------------------------------------------------------

__version__ = "1.0.a"
# Year, month, day
__last_update_date__ = "2020-04-13"

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


import os
from re import search as re_search

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
    print("\nScript 'mean_qual.py' calculates mean quality of reads in fastq file(s).\n")
    print("Version {}; {} edition.".format(__version__, __last_update_date__))
    print("\nUsage:")
    print("  python3 mean_qual.py first.fastq second.fastq.gz third.fq.gz")
    print("Following command will process all *.fastq(.gz) and *.fq(.gz) files in the working directory:")
    print("  python3 mean_qual.py")
    print("\nOptions:")
    print("  -h (--help): print help message.")
    print("  -v (--version): print version.")
    print("  -o (--outfile): output file. Default: './mean_qual_result.tsv'.")
    print("  -p (--phred-offset): Phred offset (33 of 64). Default: 33.")
    platf_depend_exit()
# end if

# Firstly check for information-providing flags

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print_help()
    platf_depend_exit()
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit()
# end if


is_fastq = lambda f: False if re_search(r".*\.f(ast)?q(\.gz)?$", f) is None else True


# Handle command line arguments
import getopt

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "p:o:",
        ["phred-offset=", "outfile="])
except getopt.GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

fpaths = args

phred_offset = 33
outfpath = "mean_qual_result.tsv"
write_head = True

for opt, arg in opts:

    if opt in ("-p", "--phred-offset"):

        try:
            phred_offset = int(arg)
            if phred_offset != 33 and phred_offset != 64:
                raise ValueError
            # end if
        except ValueError:
            print("Invalid Phred offset specified: {}".format(arg))
            print("33 and 64 are available")
            platf_depend_exit(1)
        # end try

    elif opt in ("-o", "--outfile"):
        outfpath = arg
    # end if
# end for


# If no input files are specified -- process all fastq files in the working directory
if len(fpaths) == 0:

    fpaths = tuple( filter(is_fastq, os.listdir('.')) )

    if len(fpaths) == 0:
        print("There are no *.fastq(.gz) or *.fq(.gz) files in the working directory.")
        platf_depend_exit(1)
    # end if

# If we process files in working directory,
#   their exstance and extention are already checked
else:

    # Check if all files exist
    if len(tuple(filter(os.path.exists, fpaths))) != len(fpaths):
        print("Following files do not exist:")
        for fpath in filter(lambda f: not os.path.exists(f), fpaths):
            print(fpath)
        # end for
        platf_depend_exit(1)
    # end if

    # Check if all files specified in command line are fastq files (just check extention)
    if len(tuple(filter(is_fastq, fpaths))) != len(fpaths):
        print("Following files are probably not FASTQ files:")
        for fpath in filter(lambda f: not os.path.exists(f), fpaths):
            print(fpath)
        # end for
        print("If they are actually are -- change the extention to .fastq(.gz) or .fq(.gz), please.")
        platf_depend_exit(1)
    # end if
# end if


# If output file wxists -- ask for ovewriting permissins
if os.path.exists(outfpath):

    error = True
    while error:

        reply = input("File '{}' exists. Overwrite? [y/n]:".format(outfpath))

        if reply.lower() == 'y' or reply == '':
            print("Data in file '{}' will be overwritten.\n".format(outfpath))
            error = False

        elif reply.lower() == 'n':
            print("Appending current result to '{}'.\n".format(outfpath))
            error = False
            write_head = False
        else:
            print("Invalid reply: '{}'".format(reply))
        # end if
    # end while
# end if

# Write head of the table and (maybe) overwrite old data
if write_head:
    with open(outfpath, 'w') as outfile:
        outfile.write('\t'.join( ("File",
            "Number of reads",
            "Mean quality (Q)",
            "Accuracy (%)") ) +
        '\n')
    # end with
# end if


from math import log
from gzip import open as open_as_gzip

# Function for getting Q value from Phred character:
substr_phred = lambda q_symb: ord(q_symb) - phred_offset
# List of propabilities corresponding to indices (index is Q, is the propability):
q2p_map = [10 ** (-q/10) for q in range(128)] # 127 -- max value of a signed byte
# Function for accessing propabilities by Q:
qual2prop = lambda q: q2p_map[q]
# Function for accessing Q by propability:
prop2qual = lambda p: round(-10 * log(p, 10), 2)

props = list()


for fpath in fpaths:

    print("Processing file '{}'...".format(os.path.basename(fpath)))

    if fpath.endswith(".gz"):
        open_func = open_as_gzip
        fmt_func = lambda l: l.decode("utf-8").strip()
    else:
        open_func = open
        fmt_func = lambda l: l.strip()
    # end if

    read_counter = 0

    with open_func(fpath) as fq_file:

        i = 0

        for line in fq_file:

            if i == 3:

                p = map(substr_phred, fmt_func(line))
                p = map(qual2prop, p)
                props.extend(tuple(p))

                read_counter += 1
                i = 0

            else:
                i += 1
            # end if
        # end for
    # end with

    lenpr = len(props)
    sum_props = sum(props)

    props = list()

    avg_prop = sum_props / lenpr
    accuracy = round((lenpr - sum_props) / lenpr * 100, 2)
    mean_q = prop2qual(avg_prop)

    print("\n{:,} reads".format(read_counter).replace(',', ' '))
    print("Mean quality is {}".format(mean_q))
    print("I.e. accuracy is {}%".format(accuracy))
    print('=' * 20 + '\n')

    with open(outfpath, 'a') as outfile:
        outfile.write( '\t'.join( (os.path.basename(fpath),
            str(read_counter),
            str(mean_q),
            str(accuracy)) ) +
        '\n' )
    # end with
# end for

print("Completed!")
platf_depend_exit()