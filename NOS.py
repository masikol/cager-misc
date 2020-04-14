#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# --------------------------------
# Script counts non-overalapping occurences of query sequence
#   (and it's reverse complement "comrade")
#   in all '*.fasta(.gz)' and '*.fa(.gz)' file(s).
# NOS means "non-overlapping subsequences".
# --------------------------------

import sys

if sys.version_info.major < 3:
    print( "\nYour python interpreter version is " + "%d.%d" % (sys.version_info.major,
        sys.version_info.minor) )
    print("   Please, use Python 3!\a")
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    raw_input("Press ENTER to exit:")
    sys.exit(1)
# end if

__version__ = "1.1.a"
__last_update_date__ = "2020-04-14"

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
    print("NOS; Version {}; {} edition;\n".format(__version__, __last_update_date__))
    print("""Script counts non-overalapping occurences of query sequence
  (and it's reverse complement "comrade") in '*.fasta(.gz)' and '*.fa(.gz)' file(s).""")
    print('NOS means "non-overlapping subsequences"\n')
    print("""Options:
  -h (--help) --- print help message;
  -v (--version) --- print version;
  -q (--query) --- single-sequence fasta file containing query sequence;\n""")

    print("Fasta files, in which script will search for query sequence,")
    print("  should be specified in command line without option.\n")
    print("Example:")
    print("  ./NOS.py -q query.fasta my_favourite_genome.fasta.gz another_genome.fa")
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

import os
from re import search as re_search
from getopt import gnu_getopt, GetoptError

try:
    opts, args = gnu_getopt(sys.argv[1:], "hvq:", ["help", "version", "query="])
except GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

is_fasta = lambda file: False if re_search(r"\.fa(sta)?(\.gz)?$", file) is None else True
query_fpath = None

# Get path to query file from argv
for opt, arg in opts:

    if opt in ("-q", "--query"):

        if not os.path.exists(arg):
            print("File '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        # end if

        if not is_fasta(arg):
            print("File '{}' is probably not fasta file.".format(arg))
            print("The script infers whether it's fasta file by it's extention:")
            print("'*.fasta(.gz)' or '*.fa(.gz)'")
            platf_depend_exit(1)
        # end if
        query_fpath = os.path.abspath(arg)
# end for

# Get paths to subject file(s) from argv
subject_list = args

for arg in subject_list:

    if not os.path.exists(arg):
        print("File '{}' does not exist!".format(arg))
        platf_depend_exit(1)
    # end if

    if not is_fasta(arg):
        print("File '{}' is probably not fasta file.".format(arg))
        print("The script infers whether it's fasta file by it's extention:")
        print("'*.fasta(.gz)' or '*.fa(.gz)'")
        platf_depend_exit(1)
    # end if
# end for

# If not subject file provded -- find all fasta files in working directory
if len(subject_list) == 0:
    subject_list = list(filter(is_fasta, os.listdir('.')))  # we need only fa(sta)
# end if

if len(subject_list) == 0:
    print("\nNo input data found.\n")
    print_help()
    platf_depend_exit(1)
# end if

print("NOS; Version {}; {} edition;\n".format(__version__, __last_update_date__))

# Functions for proper processing of .gz files
from gzip import open as open_as_gzip
is_gzipped = lambda f: True if f.endswith(".gz") else False

OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode("utf-8").strip()  # format gzipped line
)

# Get query
query = ""
if query_fpath is None:
    # Ask to enter query from keyboard if no query was specified in argv
    query = input("Enter query sequence (or merely press ENTER to exit):")
    if query == "":
        print("Exiting...")
        sys.exit(0)
    # end if
else:
    # Read the query
    how_to_open = OPEN_FUNCS[ is_gzipped(query_fpath) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(query_fpath) ]

    with how_to_open(query_fpath) as query_file:
        line = fmt_func(query_file.readline())
        if line == "":
            print("File '{}' is empty!".format(query_fpath))
            platf_depend_exit(1)
        # end if
        if not line.startswith('>'):
            print("File '{}' is not a valid FASTA file!".format(query_fpath))
            platf_depend_exit(1)
        # end if
        line = fmt_func(query_file.readline())
        while line:
            query += line
            line = fmt_func(query_file.readline())
        # end while
    # end with
# end if

query = query.upper()

# Configure reverse-complement query
COMPL_DICT = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'U': 'A', 'N': 'N'
}

rc_query = "".join(map(lambda n: COMPL_DICT[n], reversed(query)))

del COMPL_DICT # we do not need it any more

# Proceed
print() # just blank line

with open("NOS_result.txt", 'w') as outfile:

    outfile.write("Query sequence:\n{}\n\n".format(query))

    for fpath in subject_list:

        print("Processing '{}'...".format(fpath))
        outfile.write("=== File: '{}' ===\n\n".format(os.path.basename(fpath)))

        how_to_open = OPEN_FUNCS[ is_gzipped(fpath) ]
        fmt_func = FORMATTING_FUNCS[ is_gzipped(fpath) ]

        with how_to_open(fpath) as infile:

            line = fmt_func(infile.readline()) # read header
            seq_id = line
            line = fmt_func(infile.readline()) # read 1-st sequence line
            seq = ""
            while line:

                seq += line
                line = fmt_func(infile.readline())

                if line.startswith('>') or not line:
                    # Here we have end of current sequence
                    # Time to count occurences!

                    query_count = seq.count(query)
                    rc_query_count = seq.count(rc_query)
                    sum_count = query_count + rc_query_count

                    # Format values -- separate digit triades with spaces
                    query_count = "{:,}".format(query_count).replace(',', ' ')
                    rc_query_count = "{:,}".format(rc_query_count).replace(',', ' ')
                    sum_count = "{:,}".format(sum_count).replace(',', ' ')

                    # Write result
                    outfile.write("{}:\n".format(seq_id))
                    outfile.write("Forward query:\t{} occurences.\n".format(query_count))
                    outfile.write("Reverse query:\t{} occurences.\n".format(rc_query_count))
                    outfile.write("Totally:\t{} occurences.\n\n".format(sum_count))

                    # Reset variables
                    seq_id = line
                    seq = ""
            # end while
        # end with
    # end for
# end with

print("\nResults are in './NOS_result.txt'")
platf_depend_exit(0)