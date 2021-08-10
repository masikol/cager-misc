#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

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

__version__ = "1.2.a"
__last_update_date__ = "2021-08-10"


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
    print("""\n  |=== most-freq-subseq.py ===|\n
This script finds N most frequently occuring subsequences
   of given length for each sequence in FASTA file.""")
    print("Format of input: `*.fasta(.gz)` or `*.fa(.gz)` file and length of query subsequence.")
    print("Version: {}; {} edition.\n".format(__version__, __last_update_date__))
    print("./most-freq-subseq.py [-s <FASTA_FILE>] [-l <LENGTH>] [-h|--help]\n")
    print("Options:")
    print("  -h (--help) --- print help message;")
    print("  -v (--version) --- print versson;")
    print("  -l (--query-length) --- length of subsequence to search;")
    print("""  -s (--subject) --- FASTA file that contains sequence(s) to count
    the most frequently occuring subsequence in;""")
    print("""  -N (--num-top-occur) -- number of top the highest frequencies to show in output.
    Default value is 1 (i.e. show only the most frequently occuring subsequences);""")
    print("""  --both-strands -- count subsequences on both strands
    (default behaviour is to count only subsequences on positive strand).""")
    print("Example:")
    print("  ./most-freq-subseq.py -l 16 -s my_favourite_genome.fasta.gz -N 3")
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
from getopt import gnu_getopt, GetoptError

try:
    opts, args = gnu_getopt(
        sys.argv[1:],
        "hl:s:N:",
        [
            "help",
            "query-length",
            "subject",
            "num-top-occur",
            "both-strands"
        ]
    )
except GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

is_fasta = lambda file: not re.search(r"\.fa(sta)?(\.gz)?$", file) is None
query_len = None
sbjct_lst = list()
num_occs_print = 1
both_strands = False

for opt, arg in opts:

    if opt in ("-l", "--query-length"):

        try:
            query_len = int(arg)
            if query_len <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\aLength of query sequence must be positive integer number!")
            platf_depend_exit(1)
        # end try

    elif opt in ("-s", "--subject"):

        if not os.path.exists(arg):
            print("\aFile `{}` does not exist!".format(os.path.abspath(arg)))
            platf_depend_exit(1)
        # end if
        if not is_fasta(arg):
            print("\aFile `{}` is not a FASTA file!".format(os.path.abspath(arg)))
            platf_depend_exit(1)
        # end if
        sbjct_lst.append( os.path.abspath(arg) )

    elif opt in ("-N", "--num-top-occur"):
        try:
            num_occs_print = int(arg)
            if num_occs_print <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("""\a\nNumber of top-most-frequently-occuring subsequences
  must be positive integer number.""")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt == "--both-strands":
        both_strands = True
    # end if
# end for

if len(sbjct_lst) == 0:
    sbjct_lst = list(filter(is_fasta, os.listdir('.')))
# end if
if len(sbjct_lst) == 0:
    print("No input file found!")
    print("Please, run this script with -h flag to see help.")
    platf_depend_exit(1)
else:
    print("\nFollowing files will be processed:")
    for i, fpath in enumerate(sbjct_lst):
        print(" {}. `{}`".format(i+1, fpath))
    print('-' * 15 + '\n')
# end if


if query_len is None:
    error = True
    while error:

        query_len = input("Please, enter length of query sequence:>>")
        try:
            query_len = int(query_len)
            if query_len <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\aLength of query sequence must be positive integer number!\n")
        else:
            error = False
        # end try
    # end while
    print()
# end if


import gzip
from functools import partial
is_gzipped = lambda f: f.endswith(".gz")

OPEN_FUNCS = (
    partial(open,      mode='rt'),
    partial(gzip.open, mode='rt')
)


# Define the samt function differently
#  depending on flag `--both-strands`.
# Let's do it in order not to check this flag at every fricking iteration
if both_strands:

    def count_occurences(sequence, query_len):
        # Function counts occurences only on both strands

        # In accordance with IUPAC:
        # https://www.bioinformatics.org/sms/iupac.html
        complement_dict = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'U': 'A', 'N': 'N'
        }

        subseq_dict = dict()

        rc = lambda seq: "".join((complement_dict[nucl] for nucl in reversed(seq)))

        stop = len(sequence) - query_len

        for left in range(stop+1):

            curr_seq = sequence[left: left + query_len]

            for subseq in (curr_seq, rc(curr_seq)):

                try:
                    subseq_dict[subseq] += 1
                except KeyError:
                    subseq_dict[subseq] = 1
                # end if
            # end for
        # end while

        return subseq_dict
    # end def count_occurences

else:

    def count_occurences(sequence, query_len):
        # Function counts occurences only on positive strand

        subseq_dict = dict()

        stop = len(sequence) - query_len

        for left in range(stop+1):
            curr_seq = sequence[left: left + query_len]
            try:
                subseq_dict[curr_seq] += 1
            except KeyError:
                subseq_dict[curr_seq] = 1
            # end if
        # end while

        return subseq_dict
    # end def count_occurences
# end if


def report_N_most_freq_occs(sequence, query_len, num_occs_print, outfile):

    sequence = sequence.upper()

    # `subseq_dict` (dictionary of subsequences) has the following format:
    #   {<SUBSEQUENCE>: <OCCURRENCE>},
    # where: <SUBSEQUENCE> if of str type,
    #        <OCCURRENCE> is of int type
    subseq_dict = count_occurences(sequence, query_len)

    # Select top-'num_occs_print' entries
    subseq_list = sorted(subseq_dict.items(), key=lambda x: -x[1])
    del subseq_dict

    prev_score = subseq_list[0][1]
    N = 0

    for subseq, score in subseq_list:
        if score != prev_score:
            prev_score = score
            N += 1
        # end if
        if N == num_occs_print:
            break
        # end if
        printl("{} - {} occurences.".format(subseq, score))
    # end for
# end def report_N_most_freq_occs


print(" - Length of subsequence to search: {};".format(query_len))
print(" - Number of top-the-highest frequencies to show in output: {};".format(num_occs_print))
if both_strands:
    print(" - The script will count subsequences on both strands of input sequences.")
else:
    print(" - The script will count subsequences only on positive strand of input sequences.")
# end if
print()


with open("most_freq_subseq_result.txt", 'w') as outfile:

    def printl(text=""):
        print(text)
        outfile.write(text + '\n')
    # end def printl

    for fpath in sbjct_lst:

        printl("File `{}`:".format(fpath))
        how_to_open = OPEN_FUNCS[ is_gzipped(fpath) ]

        with how_to_open(fpath) as infile:

            sequence = ""
            line = infile.readline().strip()
            seq_id = line

            if not line.startswith('>'):
                printl("Input file is not a valid FASTA file!")
                printl("File: '{}'".format(fpath))
                platf_depend_exit(1)
            else:
                printl("Sequence {}:".format(seq_id[1:]))
            # end if

            line = infile.readline().strip()

            while line != "":

                sequence += line
                line = infile.readline().strip()

                if line.startswith('>') or line == "":

                    if query_len > len(sequence):
                        printl(
                            "Query length ({} bp) is greater then length of the subject sequence ({} bp)!" \
                                .format(query_len, len(sequence))
                        )
                    else:
                        report_N_most_freq_occs(sequence, query_len, num_occs_print, outfile)
                    # end if

                    if line.startswith('>'):
                        sequence = ""
                        seq_id = line
                        line = infile.readline().strip()
                        printl("Sequence {}:".format(seq_id[1:]))
                    # end if
                # end if
            # end while
        # end with
    # end for
# end with

platf_depend_exit(0)
