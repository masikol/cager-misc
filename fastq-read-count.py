#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# ----------------------------------------------------------------------------------
# Script counts total amount of reads in all fastq files
#-----------------------------------------------------------------------------------

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
    print("\nScript 'fastq-read-count.py' counts amount of reads in fastq file(s).\n")
    print("Version {}; {} edition.".format(__version__, __last_update_date__))
    print("\nUsage:")
    print("  python3 fastq-read-count.py first.fastq second.fastq.gz third.fq.gz")
    print("Following command will process all *.fastq(.gz) and *.fq(.gz) files in the working directory:")
    print("  python3 fastq-read-count.py")
    print("\nOptions:")
    print("  -h (--help): print help message.")
    print("  -v (--version): print version.")
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

is_fastq = lambda f: False if re_search(r".*\.f(ast)?q(\.gz)?$", f) is None else True

fpaths = list()

valid_options = ("-h", "--help", "-v", "--version")

for arg in sys.argv[1:]:

    if not arg in valid_options:

        if not os.path.exists(arg):
            print("File '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        # end if

        if not is_fastq(arg):
            print("File '{}' does not look like a fastq file".format(arg))
            print("Script understands only '*.fastq(.gz)'' and '*.fq(.gz)' extentions")
            platf_depend_exit(1)
        # end if

        fpaths.append(arg)
# end for

del valid_options

# If no files are specified in command line
if len(fpaths) == 0:
    # Get all fastq files in working dir
    fpaths = tuple(filter(is_fastq, os.listdir('.')))
# end if

if len(fpaths) == 0:
    print_help()
    platf_depend_exit(1)
# end if

print("fastq-read-count. Version {}; {} edition.\n".format(__version__, __last_update_date__))
print("Following files are found and will be processed:")
for i, f in enumerate(fpaths):
    print("{}. '{}'".format(i+1, os.path.abspath(f)))
# end for
print('-' * 45 + '\n')


import gzip

LINES_IN_READ = 4 # 4 lines per fastq record

total_read_num = 0


with open("fastq-read-count_result.tsv", 'w') as outfile:

    outfile.write('\t'.join(("File", "Number of reads")) + '\n')

    for i, fpath in enumerate(fpaths):

        read_counter = 0
        # Select open function
        open_func = gzip.open if fpath.endswith(".gz") else open

        with open_func(fpath) as infile:
            read_num = int(sum(1 for line in infile)) // LINES_IN_READ # count lines
        # end with

        total_read_num += read_num

        # Configure number with space-separated digit triades:
        str_read_num = "{:,}".format(read_num).replace(',', ' ')

        print("{}. '{}' - {} reads".format(i + 1, os.path.basename(fpath), str_read_num))
        outfile.write('\t'.join((os.path.basename(fpath), str(read_num))) + '\n')
    # end for

    # Write total number of reads
    print("Total: {:,} reads.".format(total_read_num).replace(',', ' '))
    outfile.write('\t'.join(("Total:", str(total_read_num))) + '\n')
# end with

print("Completed!")
platf_depend_exit()