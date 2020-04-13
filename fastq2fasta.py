#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# ------------------------------------------------------------
# Script converts fastq files to fasta format
#-------------------------------------------------------------

__version__ = "1.0.a"
__last_update_date__ = "2020-04-13"

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
    print("\nScript 'fastq2fasta.py' converts fastq files to fasta format.\n")
    print("Version {}; {} edition.".format(__version__, __last_update_date__))
    print("\nUsage:")
    print("  python3 fastq2fasta.py first.fastq second.fastq.gz third.fq.gz")
    print("Following command will process all *.fasta(.gz) and *.fa(.gz) files in the working directory:")
    print("  python3 fastq2fasta.py")
    print("\nOptions:")
    print("  -h (--help): print help message.")
    print("  -v (--version): print version.")
    platf_depend_exit()
# end if


# First check for information-providing flags

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print_help()
    platf_depend_exit()
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit()
# end if

fpaths = list()

from re import search as re_search

is_fastq = lambda f: False if re_search(r".*\.f(ast)?q(\.gz)?$", f) is None else True

# Add command line arguments to list of fpaths and check them
valid_options = ("-h", "--help", "-v", "--version")

for arg in sys.argv[1:]:

    if not arg in valid_options:

        if not os.path.exists(arg):
            print("File '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        # end if

        if not is_fastq(arg):
            print("File '{}' does not like an fasta file".format(arg))
            print("Scrtip understands only '*.fasta(.gz)'' and '*.fa(.gz)' extentions")
            platf_depend_exit(1)
        # end if

        fpaths.append(arg)
# end for

del valid_options

# If no input files are specified -- process all fastq files in the working directory
if len(fpaths) == 0:

    fpaths = tuple( filter(is_fastq, os.listdir('.')) )

    if len(fpaths) == 0:
        print("There are no *.fastq(.gz) or *.fq(.gz) files in the working directory.")
        platf_depend_exit(1)
    # end if
# end if

if len(fpaths) == 0:
    print_help()
# end if

print("fastq2fasta. Version {}; {} edition\n".format(__version__, __last_update_date__))

print('Following files are found and will be processed:')
for f in fpaths:
    print(os.path.abspath(f))
print('-' * 25 + '\n')

import gzip

LINES_IN_READ = 4 # 4 lines per recoed in fastq format
is_gzipped = lambda file: True if file.endswith(".gz") else False

# Start converting files:
for i, fpath in enumerate(fpaths):

    read_counter = 0

    if is_gzipped(fpath):
        # For writing bytes to file(s)
        open_func = gzip.open
        gt_chr = b'>'
        out_mode = "wb"
    else:
        # For writing strings to file(s)
        open_func = open
        gt_chr = '>'
        out_mode = "w"
    # end if

    outfpath = re_search(r"(.+)\.f(ast)?q(.gz)?$", fpath).group(1) + ".fasta"

    with open_func(fpath) as infile, open(outfpath, out_mode) as outfile:

        read_counter = 0
        # k == 1: Sequence name; k == 2: Sequence itself;
        # k == 3: Comment line; k == 4: Quality line;
        k = 1

        for line in infile:

            if k == 1: # write sequence name
                line = gt_chr + line[1:]
                outfile.write(line)

            elif k == 2: # write sequence line
                outfile.write(line)

            elif k == 4: # reset counter
                k = 0
            # end if

            k += 1
            read_counter += 1
        # end for
    # end with

    print("{}. '{}' ({} reads) --> fasta".format(i+1, os.path.basename(fpath), read_counter // LINES_IN_READ))
# end for


print("\nCompleted!")
platf_depend_exit()