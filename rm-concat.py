#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.0.a"
# Year, month, day
__last_update_date__ = "2020-11-04"

# |===== Check python interpreter version. =====|

import sys

if sys.version_info.major < 3:
    print( '\nYour python interpreter version is ' + '%d.%d' % (sys.version_info.major,
        sys.version_info.minor) )
    print('   Please, use Python 3.\a')
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith('win'):
        raw_input('Press ENTER to exit:')
    # end if
    sys.exit(1)
else:
    print('\nPython {}.{}.{}\n'.format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
# end if


import os
import re
import getopt


import gzip
# For opening plain text and gzipped files
OPEN_FUNCS = (open, gzip.open)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode('utf-8').strip()  # format gzipped line
)


from time import time, strftime, gmtime
START_TIME = time() # consider time of importing as start time
def getwt():
    """
    Function (get work time) returns time HH:MM:SS that has passed from start_time.
    """
    return strftime('%H:%M:%S', gmtime( time() - START_TIME))
# end def getwt


def platf_depend_exit(exit_code):
    """
    Function asks to press ENTER press on Windows
        and exits after that.

    :type exit_code: int;
    """
    if sys.platform.startswith('win'):
        input('Press ENTER to exit:')
    # end if
    sys.exit(exit_code)
# end def platf_depend_exit


def handle_cl_args():

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvp:', ['help', 'version' 'primers='])
    except getopt.GetoptError as opt_err:
        print( str(opt_err) )
        print('See help (`-h` option)')
        platf_depend_exit(2)
    # end try

    # Check and add fastq files:
    fq_fpaths = list()
    is_fastq = lambda f: not re.match(r'.*\.f(ast)?q(\.gz)?$', f) is None

    for arg in args:
        if not os.path.exists(arg):
            print('Error: file `{}` does not exist!'.format(arg))
            platf_depend_exit(1)
        elif not is_fastq(arg):
            print('Error: file `{}` does not look like a fastq file!'.format(arg))
            print('The script detects fastq files by extention and can process gzipped files.')
            print('Permitted extentions: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`')
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
                print('Error: file `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            elif not is_fasta(arg):
                print('Error: file `{}` does not look like a fasta file!'.format(arg))
                print('The script detects fastq files by extention and can process gzipped files.')
                print('Permitted extentions: `.fasta`, `.fa`, `.fasta.gz, `.fa.gz`')
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

    rc_dict = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G',
        'N': 'N'
    }
    single_nucl_rc = lambda nucl: rc_dict[nucl]
    rc = lambda seq: "".join(map(single_nucl_rc, seq[::-1]))

    # Retrieve primers
    primers_seqs = list()

    file_type = int(primers_fpath.endswith('.gz'))
    open_func = OPEN_FUNCS[file_type]
    fmt_func = FORMATTING_FUNCS[file_type]

    is_seq = lambda l: not l.startswith('>')
    with open_func(primers_fpath) as primers_file:
        primers_seqs = tuple(map(str.upper,                # to upper case
            map(str.strip,                                 # srtip
                filter(is_seq, primers_file.readlines())   # extract only sequences
               )
            )
        )
    # end with

    primers_seqs_and_rc = list()
    for seq in primers_seqs:
        try:
            primers_seqs_and_rc.append( (seq, rc(seq)) )
        except KeyError:
            print('Error: the script cannot handle nucleotide characters except following: ATGCN')
            print('Erroneous primer sequence: `{}`'.format(seq))
            platf_depend_exit(1)
        # end try
    # end for

    return fq_fpaths, primers_seqs_and_rc
# end def handle_cl_args


def make_outfpaths(fq_fpath):

    ok_out_fpath = os.path.join(
        os.path.dirname(fq_fpath),
        'ok-concat_{}'.format(os.path.basename(fq_fpath).replace('.gz', ''))
    )

    trash_out_fpath = os.path.join(
        os.path.dirname(fq_fpath),
        'trash-concat_{}'.format(os.path.basename(fq_fpath).replace('.gz', ''))
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
                lines.clear()
        # end for
    # end with
# end def fastq_records


def write_fastq_record(fq_record, outfile):
    outfile.write('{}\n{}\n{}\n{}\n'.format(fq_record['seq_id'],
            fq_record['seq'], fq_record['cmnt'], fq_record['qual'])
    )
# end def write_fastq_record


def concat_clean(fq_fpath, primers_seqs):

    ok_out_fpath, trash_out_fpath = make_outfpaths(fq_fpath)
    ok_count, trash_count = 0, 0

    nreads = sum(1 for _ in OPEN_FUNCS[int(fq_fpath.endswith('.gz'))](fq_fpath)) // 4
    bar_len = int(os.get_terminal_size().columns * 0.50) - 1 # minus one for >
    next_print_num = int(nreads * 0.05)
    inc_num = next_print_num

    sys.stdout.write('[{}] 0/{}'.format(' '*bar_len, nreads))
    sys.stdout.flush()

    with open(ok_out_fpath, 'w') as ok_out_file, open(trash_out_fpath, 'w') as trash_out_file:

        for i, fq_record in enumerate(fastq_records(fq_fpath)):

            read_seq = fq_record['seq'].upper()
            is_concat = False

            for seqs in primers_seqs:
                if read_seq.count(seqs[0]) + read_seq.count(seqs[1]) > 1:
                    is_concat = True
                    break
                # end if
            # end for

            if is_concat:
                write_fastq_record(fq_record, trash_out_file)
                trash_count += 1
            else:
                write_fastq_record(fq_record, ok_out_file)
                ok_count += 1
            # end if

            if i > next_print_num:
                done_ratio = i / nreads
                sys.stdout.write('\r[{}>{}] {}/{}'.format('='*int(bar_len*done_ratio),
                    ' '*int(bar_len*(1-done_ratio)), i, nreads) )
                sys.stdout.flush()
                next_print_num += inc_num
            # end if
        # end for

        sys.stdout.write('\r[{}] {}/{}\n'.format('='*bar_len, nreads, nreads))
        sys.stdout.flush()
    # end with

    print('Result files:')
    for f in (ok_out_fpath, trash_out_fpath):
        print(' `{}`'.format(f))
    # end for

    return ok_count, trash_count
# end def concat_clean


def main():
    fq_fpaths, primers_seqs = handle_cl_args()

    ok_count, trash_count = 0, 0
    print('{} - Start.'.format(getwt()))

    for fq_fpath in fq_fpaths:
        print('{} - Processing file `{}`.'.format(getwt(), fq_fpath))
        ok_count_tmp, trash_count_tmp = concat_clean(fq_fpath, primers_seqs)
        ok_count += ok_count_tmp
        trash_count += trash_count_tmp
        print('{} - File `{}` is procesed.'.format(getwt(), fq_fpath))
        print('-'*10)
    # end for

    print('Results: {} reads keeped, {} reads discarded.'.format('{:,}'.format(ok_count).replace(',', ' '),
                                                                 '{:,}'.format(trash_count).replace(',', ' '))
    )
# end def main


if __name__ == '__main__':
    main()
# end if