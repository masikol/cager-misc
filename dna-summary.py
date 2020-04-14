#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

__version__ = "1.0.a"
__last_update_date__ = "2020-04-13"

# ---------------------------------------------------------------------------
# Script makes brief summary (length, coverage, GC-content)
#   of .dna files located in directory './contigs'
#   and saves it in file 'dna-summary.txt' in tab-separated format.
#
# Coverage can be extracted if file is named as SPAdes contig:
#   NODE_1_length_61704_cov_114.517.dna
#
# Format of output file:
# Ordinal number  <tab>  File <tab> Length (b.p.) <tab> Coverage <tab> GC-content (%)
# 1	NODE_1_length_61704_cov_114.517.dna	61704	114,517 	44,54
#
# ---------------------------------------------------------------------------


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
    print("\n    |=== dna-summary ===|")
    print("Version {}. {} edition.\n".format(__version__, __last_update_date__))

    print("""Script makes brief summary (length, coverage, GC-content)
  of .dna files located in directory './contigs'
  and saves it in file 'dna-summary.txt' in tab-separated format.""")

    print("""Coverage can be extracted if file is named as SPAdes contig:
  NODE_1_length_61704_cov_114.517.dna""")
    print("""Format of output file:
  Ordinal number  <tab>  File <tab> Length (b.p.) <tab> Coverage <tab> GC-content (%)
  1   NODE_1_length_61704_cov_114.517.dna 61704   114,517     44,54""")

    print("\nUsage: just run it in directory, where folder 'contigs' is located")
    print("  python3 dna-summary.py")
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

import os
import re

# Specially for evil russian KGB
import locale
if locale.getdefaultlocale()[0].startswith('ru'):
    dec_sep = ','
else:
    dec_sep = '.'
# end if

target_dir = "contigs"

# Check if there is 'contigs/' directory in working dir
if not os.path.exists(target_dir):
    print('\n"contigs/" directory is not found!\n')
    print_help()
    platf_depend_exit(1)
# end if

print("dna-summary; Version {}; {} edition;".format(__version__, __last_update_date__))
print('"contigs/" directory is found and will be processed.\n')

# We will process only those .dna files, which are SPAdes contigs
is_target_file = lambda f: not re.match(r"NODE_[0-9]+_length_[0-9]+_cov_[0-9\.\,]+\.dna", f) is None


fpaths = os.listdir(target_dir) # get all files in 'contigs' dir
fpaths = list(filter(is_target_file, fpaths)) # leave only 'NODE' files

if len(fpaths) == 0:
    print("There is no 'NODE...' files in 'contigs' directory")
    print_help()
    platf_depend_exit(1)
# end if

# Sort the list by NODE ordinal number
fpaths.sort(key=lambda x: int(x.split('_')[1]))

# Variables for summary
total_length = 0
mean_coverage = 0
min_coverage = float('inf')
max_coverage = 0

# Start processing .dna files
with open('dna-summary.txt', 'w') as outfile:

    outfile.write('\t'.join( ("#",
        "Sequence name",
        "Length (b.p.)",
        "Coverage",
        "GC (%)") ) +
    '\n')

    for i, fpath in enumerate(fpaths):

        length = int(fpaths[i].split('_')[3]) # get length

        # Parse .dna file (read all):
        with open(os.path.join(target_dir, fpath), 'rb') as infile:
            seq = infile.read()
        # end with

        # Trim leading 25 bytes (till the sequence starts), extract the sequence and decode it
        seq = seq[25:25+length].decode('ascii', 'ignore')

        # Count GC-content
        g_count = seq.count('G')
        c_count = seq.count('C')
        s_count = seq.count('S')
        gc_content = str(round(((g_count + c_count + s_count) / length * 100), 2))

        cov = float(fpath.split('_')[5].replace(".dna", '')) # get coverage

        cov_str = str(cov).replace('.', dec_sep) # evil russian KGB
        gc_content = gc_content.replace('.', dec_sep) # evil russian KGB

        outfile.write('\t'.join( (i+1,
            fpath.replace(".dna", ''),
            str(len(seq)),
            cov_str,
            gc_content) ) +
        '\n')

        # Collect some info
        total_length += length
        mean_coverage += cov
        max_coverage = max(max_coverage, cov)
        min_coverage = min(min_coverage, cov)
    # end for

    outfile.write("\n")
    # Add sone extra data
    for print_func in (sys.stdout.write, outfile.write):
        # Format length -- separate digit triades with spaces
        print_func("Total length:   {:,}\n".format(total_length).replace(',', ' '))
        print_func('Min Coverage:   {}\n'.format(min_coverage))
        print_func('Max Coverage:   {}\n'.format(max_coverage))
        print_func('Mean Coverage:  {}\n'.format(round(mean_coverage / len(fpaths), 3)))
        print_func('Files processed: {}\n'.format(len(fpaths)))
    # end for

# end with
print('Completed!')
platf_depend_exit()