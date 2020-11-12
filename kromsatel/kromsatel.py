#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.2.a"
# Year, month, day
__last_update_date__ = "2020-11-12"

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
import subprocess as sp

import io

try:
    import numpy as np
except ImportError:
    print('`numpy` package is not installed')
    print('Please, install it. Example command:')
    print('  pip3 install numpy')
    platf_depend_exit(1)
# end try

try:
    import pandas as pd
except ImportError:
    print('`pandas` package is not installed')
    print('Please, install it. Example command:')
    print('  pip3 install pandas')
    platf_depend_exit(1)
# end try

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
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hd:t:', ['help', 'db=', 'threads='])
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

    if '-h' in sys.argv[1:] or '--help' in sys.argv[1:] or '-help' in sys.argv[1:]:
        print('kromsatel')
        print('Version {}. {} edition.'.format(__version__, __last_update_date__))
        print("Usage:")
        print('  ./kromsatel.py <YOUR_READS> -d fragments-db/nCoV-2019_ref-fragments.fasta')
        print('For example:')
        print('  ./kromsatel.py <YOUR_READS> -d fragments-db/nCoV-2019_ref-fragments.fasta')
        platf_depend_exit(0)
    # end if

    if len(fq_fpaths) == 0:
        print("Usage:")
        print('./kromsatel.py <YOUR_READS> -d fragments-db/nCoV-2019_ref-fragments.fasta')
        platf_depend_exit(0)
    # end if

    db_fpath = None
    n_thr = 1

    # Find path to db
    for opt, arg in opts:
        if opt in ('-d', '--db'):
            if not os.path.exists( '{}.nhr'.format(arg) ):
                print('Error: database `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            else:
                db_fpath = os.path.abspath(arg)
            # end if

        elif opt in ('-t', '--threads'):
            try:
                n_thr = int(arg)
                if n_thr < 1:
                    raise ValueError
                # end if
            except ValueError:
                print('Error: number of threads must be positive integer number!')
                print(' And here is your value: `{}`'.format(arg))
                sys.exit(1)
            # end try
            if n_thr > len(os.sched_getaffinity(0)):
                print('''\nWarning! You have specified {} threads to use
      although {} are available.'''.format(n_thr, len(os.sched_getaffinity(0))))
                error = True
                while error:
                    reply = input('''\nPress ENTER to switch to {} threads,
      or enter 'c' to continue with {} threads,
      or enter 'q' to exit:>>'''.format(len(os.sched_getaffinity(0)), n_thr))
                    if reply in ('', 'c', 'q'):
                        error = False
                        if reply == '':
                            n_thr = len(os.sched_getaffinity(0))
                            print('\nNumber of threads switched to {}\n'.format(n_thr))
                        elif reply == 'c':
                            pass
                        elif reply == 'q':
                            sys.exit(0)
                        # end if
                    else:
                        print('\nInvalid reply!\n')
                    # end if
                # end while
            # end if
            # end if
        # end if
    # end for

    return fq_fpaths, db_fpath, n_thr
# end def handle_cl_args


def fastq2fasta(fq_fpath):

    read_counter = 0
    LINES_IN_READ = 4

    if fq_fpath.endswith('.gz'):
        # For writing bytes to file(s)
        open_func = gzip.open
        gt_chr = b'>'
        out_mode = "wb"
        endl = b'\n'
        space = b' '
    else:
        # For writing strings to file(s)
        open_func = open
        gt_chr = '>'
        out_mode = "w"
        endl = '\n'
        space = ' '
    # end if

    outfpath = re.search(r"(.+)\.f(ast)?q(.gz)?$", fq_fpath).group(1) + ".fasta"

    with open_func(fq_fpath) as infile, open(outfpath, out_mode) as outfile:

        read_counter = 0
        # k == 1: Sequence name; k == 2: Sequence itself;
        # k == 3: Comment line; k == 4: Quality line;
        k = 1

        for line in infile:

            if k == 1: # write sequence name
                line = gt_chr + line[1:].partition(space)[0]
                outfile.write(line + endl)

            elif k == 2: # write sequence line
                outfile.write(line)

            elif k == 4: # reset counter
                k = 0
            # end if

            k += 1
            read_counter += 1
        # end for
    # end with

    print('{} - `{}` ({} reads) --> fasta'.format(getwt(), os.path.basename(fq_fpath), read_counter // LINES_IN_READ))

    return outfpath
# end def fastq2fasta


def disco_align(query_fpath, db_fpath, n_thr):

    # Check if 'blast+' tookit is installed
    pathdirs = os.environ['PATH'].split(os.pathsep)
    utility = 'blastn'
    utility_found = False
    for directory in pathdirs:
        if os.path.exists(directory) and utility in os.listdir(directory):
            utility_found = True
            break
        # end if
    # end for
    if not utility_found:
        print('  Attention!\n`{}` from BLAST+ toolkit is not installed.'.format(utility))
        print('''If this error still occures although you have installed everything 
-- make sure that this program is added to PATH)''')
        platf_depend_exit(1)
    # end if

    outfmt = '6 qseqid sseqid sstrand length qlen slen qstart qend sstart send'

    # Configure command line
    blast_cmd = "blastn -query {} \
    -db {} \
    -task {} \
    -num_threads {} \
    -evalue 1e-3 \
    -outfmt '{}'".format(query_fpath,
                         db_fpath,
                         'dc-megablast',
                         n_thr,
                         outfmt
    )

    sys.stdout.write('{} - Aligning...'.format(getwt()))
    sys.stdout.flush()

    # Launch blastn
    pipe = sp.Popen(blast_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('\nError while aligning a sequence against local database')
        print(stdout_stderr[1].decode('utf-8'))
        platf_depend_exit(pipe.returncode)
    # end if

    sys.stdout.write('\r{} - Aligning...done\n'.format(getwt()))
    sys.stdout.flush()

    return stdout_stderr[0].decode('utf-8')
# end def disco_align


def str2df(align_result_str):

    # aln_df = pd.DataFrame([l.split('\t') for l in align_result_str.split('\n')])
    aln_df = pd.read_csv(io.StringIO(align_result_str),
        sep='\t',
        header=None,
        names= ['qseqid', 'sseqid', 'sstrand',
                'length',   'qlen',    'slen',
                'qstart',   'qend',  'sstart', 
                                       'send'],
        true_values=['plus'], false_values=['minus']
    )

    return aln_df
# end def str2df


def make_outfpath(fq_fpath):

    return os.path.join(
        os.path.dirname(fq_fpath),
        'cleaned_{}'.format(os.path.basename(fq_fpath).replace('.gz', ''))
    )
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
                    'seq_id': lines[0][1:].partition(' ')[0],
                    'seq': lines[1],
                    'cmnt': lines[2],
                    'qual': line
                }
                lines.clear()
        # end for
    # end with
# end def fastq_records


def write_fastq_record(fq_record, outfile):
    outfile.write('@{}\n{}\n{}\n{}\n'.format(fq_record['seq_id'],
            fq_record['seq'], fq_record['cmnt'], fq_record['qual'])
    )
# end def write_fastq_record


# Positions in BLAST+ are 1-based
MAX_EDGE_OFFSET = 1


def set_touch_start(row):

    if row['sstrand']:
        row['touch_start'] = row['sstart'] < MAX_EDGE_OFFSET + 2
    else:
        row['touch_start'] = row['sstart'] > row['slen'] - 2
    # end if

    return row
# end def get_touch_start


def set_touch_end(row):

    if row['sstrand']:
        row['touch_end'] = row['send'] > row['slen'] - 1 - MAX_EDGE_OFFSET
    else:
        row['touch_end'] = row['send'] < MAX_EDGE_OFFSET + 2
    # end if

    return row
# end def get_touch_start


def map_ends(curr_alns):

    if curr_alns.empty:
        return []
    # end if

    curr_alns = curr_alns.apply(set_touch_start, axis=1)
    curr_alns = curr_alns.apply(set_touch_end, axis=1)

    full_span_alns = curr_alns.query('touch_start & touch_end')
    one_side_alns = curr_alns.query('touch_start | touch_end')

    aligned_fragments = list()

    # Get major fragment first.
    # Their accessions start with 'A'
    full_span_alns = full_span_alns.sort_values(by='sseqid', ascending=True)
    # full_span_alns.sort_values(by='sseqid', ascending=True, inplace=True)

    cov_array = np.zeros(curr_alns.iloc[0,:]['qlen'], dtype=np.uint8)


    for aln_collenstion in (full_span_alns, one_side_alns):

        for i in range(aln_collenstion.shape[0]):

            aln = aln_collenstion.iloc[i, :]
            buff_cov_array = np.copy(cov_array)
            np.add.at(buff_cov_array, range(aln['qstart']-1, aln['qend']), 1)

            if not any(buff_cov_array > 1):
                np.add.at(cov_array, range(aln['qstart']-1, aln['qend']), 1)
                aligned_fragments.append( (aln['qstart']-1, aln['qend']) )
            # end if
        # end for
    # end for

    return aligned_fragments

# end def map_ends


def clean_and_shred(aln_df, fq_fpath):

    outfpath = make_outfpath(fq_fpath)

    nreads = sum(1 for _ in OPEN_FUNCS[int(fq_fpath.endswith('.gz'))](fq_fpath)) // 4
    bar_len = int(os.get_terminal_size().columns * 0.50)
    next_print_num = int(nreads * 0.01)
    inc_num = next_print_num

    print('{} - Start cleaning and shredding...'.format(getwt()))

    sys.stdout.write('{} - [{}] 0/{}'.format(getwt(), ' '*bar_len, nreads))
    sys.stdout.flush()

    with open(outfpath, 'w') as outfile:

        for i, fq_record in enumerate(fastq_records(fq_fpath)):

            # curr_alns = aln_df[aln_df['qseqid'] == fq_record['seq_id']]
            curr_alns = aln_df.query('qseqid == @fq_record["seq_id"]')
            aligned_fragments = map_ends(curr_alns)

            for fragm in aligned_fragments:

                qstart, qend = fragm[0], fragm[1]
                curr_seq_id = '{}_{}-{}'.format(fq_record['seq_id'], qstart+1, qend) # write 1-based coordinates here
                curr_seq = fq_record['seq'][qstart : qend]
                curr_qual = fq_record['qual'][qstart : qend]

                write_fastq_record({'seq_id':  curr_seq_id,
                                    'seq'   :  curr_seq,
                                    'cmnt'  :  fq_record['cmnt'],
                                    'qual'  :  curr_qual},
                                   outfile)
            # end for

            if i > next_print_num:
                 bar_len = int(os.get_terminal_size().columns * 0.50)
                 done_ratio = i / nreads
                 sys.stdout.write('\r{} - [{}>{}] {}/{} ({}%)'.format(getwt(),
                    '='*(int(bar_len*done_ratio)-1),
                     ' '*int(bar_len*(1-done_ratio)), i, nreads, int(done_ratio*100)) )
                 sys.stdout.flush()
                 next_print_num += inc_num
             # end if
        # end for
    # end with

    sys.stdout.write('\r{} - [{}] {}/{} (100%)\n'.format(getwt(),
        '='*bar_len, nreads, nreads))
    sys.stdout.flush()

    return outfpath
# end def clean_and_shred


def main():
    fq_fpaths, db_fpath, n_thr = handle_cl_args()

    print('{} - Start.'.format(getwt()))

    for fq_fpath in fq_fpaths:
        print('{} - Processing file `{}`'.format(getwt(), fq_fpath))

        query_fpath = fastq2fasta(fq_fpath)

        aln_df = str2df(disco_align(query_fpath, db_fpath, n_thr))

        outfpath = clean_and_shred(aln_df, fq_fpath)

        try:
            os.unlink(query_fpath)
        except OSError as oserr:
            print('Warning: Cannot remove temporary file `{}`'.format(query_fpath))
            print( str(oserr) )
        # end try

        print('{} - File `{}` is processed.'.format(getwt(), fq_fpath))
        print('Output file: `{}`'.format(outfpath))
        print('-------')
    # end for
# end def main


if __name__ == '__main__':
    main()
# end if