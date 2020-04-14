#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
# ---------------------------------
# Script packs all "SPAdes-like" .dna files in 'contigs' directory
#   to single multi-fasta file 'packed_dna.fasta'.
# "SPAdes-like" means that name of file if of following format:
#   NODE_1_length_61704_cov_114.517.dna
# ---------------------------------

__version__ = "1.0.a"
__last_update_date__ = "2020.04.14"


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
    print("\n    |=== packer-dna-to-fasta ===|")
    print("Version {}. {} edition.\n".format(__version__, __last_update_date__))

    print('Script packs all "SPAdes-like" .dna files in \'contigs\' directory')
    print("  to single multi-fasta file 'packed_dna.fasta'.")
    print('SPAdes-like" means that name of file if of following format:')
    print("  NODE_1_length_61704_cov_114.517.dna")

    print("\nUsage: just run it in directory, where folder 'contigs' is located")
    print("  python3 packer-dna-to-fasta.py")
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

target_dir = "contigs"

# Check if there is 'contigs/' directory in working dir
if not os.path.exists(target_dir):
    print('\n"contigs/" directory is not found!\n')
    print_help()
    platf_depend_exit(1)
# end if

print("packer-dna-to-fasta; Version {}; {} edition;".format(__version__, __last_update_date__))
print('\n"contigs/" directory is found and will be processed.\n')

# We will process only those .dna files, which are SPAdes contigs
import re
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

outfpath = "packed_dna.fasta"

# Proceed
with open(outfpath, 'w') as outfile:

    for fpath in fpaths:

        length = int(fpath.split('_')[3]) # get length

        # Parse .dna file (read all):
        with open(os.path.join(target_dir, fpath), 'rb') as infile:
            seq = infile.read()
        # end with

        # Trim leading 25 bytes (till the sequence starts), extract the sequence and decode it
        seq = seq[25:25+length].decode('ascii', 'ignore')

        header = os.path.basename(fpath).replace('.dna', '')
        outfile.write(">{}\n{}\n".format(header, seq))
        print("Sequence {} is added to '{}'".format(header, outfpath))
    # end for
# end with

# Print some summary
print('=' * 45)
print('Completed!')
print("{} sequences have been written to file '{}'".format(len(fpaths), outfpath))
platf_depend_exit(0)