#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
# ---------------------------------
# Script moves "SPAdes-like" *.dna files with coverage less than specified one
#   from 'contigs/' directory to directory 'cov_below_x/'.
# "SPAdes-like" means that name of file is of following format:
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
    print("\n    |=== seqator ===|")
    print("Version {}. {} edition.\n".format(__version__, __last_update_date__))

    print("Script moves \"SPAdes-like\" *.dna files with coverage less than specified one")
    print("  from 'contigs/' directory to directory 'cov_below_x/'.")
    print("\"SPAdes-like\" means that name of file is of following format:")
    print("  NODE_1_length_61704_cov_114.517.dna")

    print("\nUsage:")
    print("You can specify coverage threshold in command line:\n")
    print("  python3 seqator.py 12\n")
    print("Or run just this, and script will ask you to enter threshold from keyboard:")
    print("  python3 seqator.py")
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

# Get coverage value
if len(sys.argv) == 2:
    # If it is specified in argv

    try:
        cov_threshold = float(sys.argv[1].replace(',', '.')) # specially for evil russian KGB
        if cov_threshold < 0:
            raise ValueError
    except ValueError:
        print("Invalid coverage specified: '{}'".format(cov_threshold))
        print("Must be positive number!")
        platf_depend_exit(1)
    # end try

elif len(sys.argv) == 1:

    # If no coverage specified in argv, ask to enter it from keyboard
    cov_threshold = None
    while cov_threshold is None:
        cov_threshold = input("Please, enter coverage threshold:")
        cov_threshold = cov_threshold.replace(',', '.') # specially for evil russian KGB
        try:
            cov_threshold = float(cov_threshold)
            if cov_threshold < 0:
                raise ValueError
        except ValueError:
            print("Invalid coverage specified: '{}'".format(cov_threshold))
            print("Must be positive number!")
            cov_threshold = None
        # end try
    # end while

else:
    print_help()
    platf_depend_exit(1)
# end if

import os

target_dir = "contigs"

# Check if there is 'contigs/' directory in working dir
if not os.path.exists(target_dir):
    print('\n"contigs/" directory is not found!\n')
    print_help()
    platf_depend_exit(1)
# end if

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

print("packer-dna-to-fasta; Version {}; {} edition;".format(__version__, __last_update_date__))
print('\n"contigs/" directory is found and will be processed.\n')

new_path = os.path.join(target_dir, "cov_below_{}".format(cov_threshold))
if not os.path.isdir(new_path):
    os.makedirs(new_path)
# end if
print("Sequences with coverage less than {} will be moved to directory '{}'\n".format(cov_threshold, new_path))

import shutil

# Proceed
moved_counter = 0
for fpath in fpaths:

    # Retrieve coverage from file name
    cov = fpath.split('_')[5].replace(',', '.')      # Get coverage + .dna
    cov = cov.replace('.dna', '')  # Prune '.dna'

    # Compare and move if less
    if float(cov) < cov_threshold:
        print("Moving '{}'".format(fpath))
        shutil.move(os.path.join(target_dir, fpath), new_path)
        moved_counter += 1
    # end if
# end for

if len(os.listdir(new_path)) == 0:
    os.rmdir(new_path)
# end if

print("=" * 45)
print("Completed!")
print("{} files moved to '{}'".format(moved_counter, new_path))
print("{} files left in '{}'".format(len(fpaths) - moved_counter, target_dir))
platf_depend_exit(0)