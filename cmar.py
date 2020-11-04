#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.0.a"
# Year, month, day
__last_update_date__ = "2020-11-04"

# |===== Check python interpreter version. =====|

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
else:
    print("\nPython {}.{}.{}\n".format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
# end if


import os
import re
import getopt
import difflib


import gzip
# For opening plain text and gzipped files
OPEN_FUNCS = (open, gzip.open)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode("utf-8").strip()  # format gzipped line
)


from time import time, strftime, localtime, gmtime
START_TIME = time() # consider time of importing as start time
def getwt():
    """
    Function (get work time) returns time HH:MM:SS that has passed from start_time.
    """
    return strftime("%H:%M:%S", gmtime( time() - START_TIME))
# end def getwt


def platf_depend_exit(exit_code):
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


def handle_cl_args():

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvp:', ['help', 'version' 'primers='])
    except getopt.GetoptError as opt_err:
        print( str(opt_err) )
        print("See help ('-h' option)")
        platf_depend_exit(2)
    # end try

    # Check and add fastq files:
    fq_fpaths = list()
    is_fastq = lambda f: not re.match(r'.*\.f(ast)?q(\.gz)?$', f) is None

    for arg in args:
        if not os.path.exists(arg):
            print("Error: file '{}' does not exist!".format(arg))
            platf_depend_exit(1)
        elif not is_fastq(arg):
            print("Error: file '{}' does not look like a fastq file!".format(arg))
            print('The script detects fastq files by extention and can process gzipped files.')
            print("Permitted extentions: '.fastq', '.fq', '.fastq.gz, '.fq.gz'")
            platf_depend_exit(1)
        else:
            fq_fpaths.append(os.path.abspath(arg))
        # end if
    # end for

    primers_fpath = None
    is_fasta = lambda f: not re.match(r'.*\.f(ast)?a(\.gz)?$', f) is None

    # Handle file containing primers
    for opt, arg in opts:
        if opt in ('-p', '--primers'):
            if not os.path.exists(arg):
                print("Error: file '{}' does not exist!".format(arg))
                platf_depend_exit(1)
            elif not is_fasta(arg):
                print("Error: file '{}' does not look like a fasta file!".format(arg))
                print('The script detects fastq files by extention and can process gzipped files.')
                print("Permitted extentions: '.fasta', '.fa', '.fasta.gz, '.fa.gz'")
                platf_depend_exit(1)
            else:
                primers_fpath = os.path.abspath(arg)
            # end if
        # end if
    # end for

    if primers_fpath is None:
        # Try to find file with primers in the working directory
        cw_fpaths = filter(os.path.isfile, os.listdir('.'))
        primer_pattern = r'.*[Pp]rimers[^{}]*\.f(ast)?a(\.gz)?$'.format(os.sep)
        prob_primer_fpaths = tuple(filter(lambda f: re.match(primer_pattern, f), cw_fpaths))

        if len(prob_primer_fpaths) == 0:
            print('Please specify fasta file containing primer sequences with `-p` option.')
            platf_depend_exit(1)
        elif len(prob_primer_fpaths) == 1:
            primers_fpath = prob_primer_fpaths[0]
            error = True
            while error:
                print('Found fasta file containing primers: `{}`'.format(primers_fpath))
                reply = input('Use it? [Y/n]:')
                if reply == '' or reply.lower() == 'y':
                    print('Using `{}` as file containing primer sequence.'.format(primers_fpath))
                    error = False
                elif reply.lower() == 'n':
                    print('Well, see you soon!')
                    platf_depend_exit(1)
                else:
                    print('Invalid reply: `{}`\n'.format(reply))
                # end if
            # end while
        else:
            print('Found multiple files probably containing primer sequences.')
            print('Here they are:')
            for i, f in enumerate(prob_primer_fpaths):
                print(' {}. `{}`'.format(i+1, f))
            # end for
            print('Please choose one with `-p` option')
            platf_depend_exit(1)
        # end if
    # end if

    # Retrieve primers
    primers_seqs = list()

    file_type = int(primers_fpath.endswith('.gz'))
    open_func = OPEN_FUNCS[file_type]
    fmt_func = FORMATTING_FUNCS[file_type]

    is_seq = lambda l: not l.startswith('>')
    with open_func(primers_fpath) as primers_file:
        primers_seqs = tuple(map(str.strip, filter(is_seq, primers_file.readlines())))
    # end with

    return fq_fpaths, primers_seqs
# end def handle_cl_args


def make_outfpaths(fq_fpath):

    ok_out_fpath = os.path.join(
        os.path.dirname(fq_fpath),
        'ok-cmar_{}'.format(os.path.basename(fq_fpath))
    )

    trash_out_fpath = os.path.join(
        os.path.dirname(fq_fpath),
        'trash-cmar_{}'.format(os.path.basename(fq_fpath))
    )

    return ok_out_fpath, trash_out_fpath
# end def make_outfpath


def fastq_records(fq_fpath):

    file_type = int(fq_fpath.endswith('.gz'))
    open_func = OPEN_FUNCS[file_type]
    fmt_func = FORMATTING_FUNCS[file_type]

    cnt = 0
    lines = list()
    with open_func(fq_fpath) as fq_file:

        for line in fq_file:

            line = fmt_func(line)

            if line == '':
                continue
            # end if

            if cnt != 3:
                lines.append(line)
                cnt += 1
            else:
                cnt = 0
                yield {
                    'seq_id': lines[0],
                    'seq': lines[1],
                    'cmnt': lines[2],
                    'qual': line
                }
        # end for
    # end with
# end def fastq_records


def write_fastq_record(fq_record, outfile):
    outfile.write('{}\n{}\n{}\n{}\n'.format(fq_record['seq_id'],
                                            fq_record['seq'],
                                            fq_record['cmnt'],
                                            fq_record['qual'])
    )
# end def write_fastq_record


def cmar_clean(fq_fpath, primers_seqs):

    ok_out_fpath, trash_out_fpath = make_outfpaths(fq_fpath)

    matcher = difflib.SequenceMatcher(isjunk = lambda x: x == 'N')

    with open(ok_out_fpath, 'w') as ok_out_file, open(trash_out_fpath, 'w') as trash_out_file:
        for fq_record in fastq_records(fq_fpath):
            print('read')
            matcher.set_seq1(fq_record['seq'].upper())
            for prim_seq in primers_seqs:
                print('primer')
                matcher.set_seq2(prim_seq.upper())

                print(matcher.get_matching_blocks())
                input()

            if is_mar:
                write_fastq_record(fq_record, trash_out_fpath)
            else:
                write_fastq_record(ok_out_fpath, trash_out_fpath)
        # end for
    # end with
# end def cmar_clean


def main():
    fq_fpaths, primers_seqs = handle_cl_args()

    print('{} - Start processing'.format(getwt()))
    for fq_fpath in fq_fpaths:
        print('{} - Processing file `{}`.'.format(getwt(), fq_fpath))
        cmar_clean(fq_fpath, primers_seqs)
        print('{} - File `{}` is procesed.'.format(getwt(), fq_fpath))
    # end for
# end def main

if __name__ == '__main__':
    main()
# end if