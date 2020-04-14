#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# -------------------------------
# Script finds fasta record(s) in fasta file by given sequence header and prints it(them) to stdout.
# -------------------------------

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

__version__ = "1.0.a"
__last_update_date__ = "2020-04-14"


# Firstly check for information-providing flags

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print_help()
    sys.exit(0)
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    sys.exit(0)
# end if


def print_help():

    print("find-seq; Version {}; {} edition;\n".format(__version__, __last_update_date__))
    print("Script finds fasta record(s) in *.fasta(.gz) or *.fa(.gz) file by given sequence header")
    print("  and prints found record(s) to stdout.")
    print("-h (--help): print help message;\n-v (--version): print version.\n")
    print("Usage: to find all fasta records in file 'some_seqs.fasta'")
    print("  whose headers contain string 'length_114', you should run:\n")
    print("  python3 find-seq.py length_114 some_seqs.fasta\n")
# end def print_help


# Scheck if there ar eonly 2 arguments in command line
if len(sys.argv[1:]) != 2:
    print_help()
    sys.exit(1)
# end if

import os
import re

is_fasta = lambda file: False if re.search(r"\.f(ast)?a(\.gz)?$", file) is None else True

# Path to file is 2-d argument
fpath = sys.argv[2]

# Check existanse
if not os.path.exists(fpath):
    print("File '{}' does not exist!".format(fpath))
    sys.exit(1)
# end if

# Check if it is fasta file
if not is_fasta(fpath):
    print("File '{}' is probably not fasta file.".format(fpath))
    print("The script infers whether it's fasta file by it's extention:")
    print("'*.fasta(.gz)' or '*.fa(.gz)'")
    sys.exit(1)
# end if

# Configure functions and variales for proper processing of .gz files

if fpath.endswith(".gz"):

    import gzip

    open_func = gzip.open   # open function
    seq_prefix = b'>'       # fasta header indicator
    # We should decode bytes to write it to stdout
    decode_func = lambda x: x.decode("utf-8")
    # Encode query in order not to decode bytes every time
    query_seq = sys.argv[1].encode("utf-8")
else:
    open_func = open        # open function
    seq_prefix = '>'        # fasta header indicator
    # If input file is plain fasta, we do not need to decode what we want to write to stdout
    decode_func = lambda x: x
    # Do not encode query string -- just get it from argv
    query_seq = sys.argv[1]
# end if

# Proceed
with open_func(fpath) as infile:

    line = infile.readline()

    while line:

        # If we find query string in any header
        if line.startswith(seq_prefix) and query_seq in line:

            # Write this header
            sys.stdout.write(decode_func(line))
            line = infile.readline() # get next line

            # Write file's content until end of file or next header
            while line and not line.startswith(seq_prefix):
                sys.stdout.write(decode_func(line))
                line = infile.readline()
            # end while
        else:
            line = infile.readline() # get next line
        # end if
    # end for
# end with

sys.exit(0)