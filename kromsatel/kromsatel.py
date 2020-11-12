#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.2.b"
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


def platf_depend_exit(exit_code):
    # Function asks to press ENTER press on Windows
    #     and exits after that.

    # :type exit_code: int;

    if sys.platform.startswith('win'):
        input('Press ENTER to exit:')
    # end if
    sys.exit(exit_code)
# end def platf_depend_exit


# Firstly check if ve just need to print version or help message
if '-v' in sys.argv[1:] or '--version' in sys.argv[1:] or '-version' in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if

if '-h' in sys.argv[1:] or '--help' in sys.argv[1:] or '-help' in sys.argv[1:]:
    print('kromsatel')
    print('Version {}. {} edition.'.format(__version__, __last_update_date__))
    print("Usage:")
    print('  ./kromsatel.py <YOUR_READS> -d <DATABASE_WITH_FRAGMENTS>')
    print('For example:')
    print('  ./kromsatel.py corona_reads.fastq -d fragments-db/nCoV-2019_ref-fragments.fasta')
    platf_depend_exit(0)
# end if

import os
import re
import getopt
import subprocess as sp

import operator
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
    # Function (get work time) returns time HH:MM:SS that has passed from start_time.
    return strftime('%H:%M:%S', gmtime( time() - START_TIME))
# end def getwt


def handle_cl_args():
    # Function handles command line arguments.
    # Returns:
    #  1. List of paths to fastq files to be processed.
    #  2. Path to reference database with fragments.
    #  3. Number of threads to launch.

    # Get arguments
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hvd:t:', ['help', 'version',
                                                                'db=', 'threads=']
        )
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

    # Check if there are any files to process
    if len(fq_fpaths) == 0:
        print('Usage:')
        print('./kromsatel.py <YOUR_READS> -d <DATABASE_WITH_FRAGMENTS>')
        platf_depend_exit(0)
    # end if

    db_fpath = None
    n_thr = 1

    # Handle options
    for opt, arg in opts:
        # Get path to database
        if opt in ('-d', '--db'):
            if not os.path.exists( '{}.nhr'.format(arg) ):
                print('Error: database `{}` does not exist!'.format(arg))
                platf_depend_exit(1)
            else:
                db_fpath = os.path.abspath(arg)
            # end if

        # Handle number of threads
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
    # end for

    return fq_fpaths, db_fpath, n_thr
# end def handle_cl_args


def fastq2fasta(fq_fpath):
    # Function converts fastq file to fasta format.
    #
    # :param fq_fpath: path to fastq file to be converted;
    # :type fq_fpath: str;
    #
    # Returns path to result fasta file.

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
    # Function alignes reads against fragments using Discontiguous Megablast.
    #
    # :param query_fpath: path to fasta query file;
    # :type query_fpath: str;
    # :param db_fpath: path to database containing target fragments;
    # :type db_fpath: str;
    # :param n_thr: number of threads to launch;
    # :type n_thr: int;
    #
    # Returns tabular string if format "outfmt 6" returned by blastn.

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
    # Function converts tabular string data to pandas.DataFrame.
    #
    #:param align_result_str: tabular string to be converted to pandas.DataFrame;
    #:type align_result_str: str;
    #
    # Returns pandas.DataFrame contatining the same information as the string passed.

    # aln_df = pd.DataFrame([l.split('\t') for l in align_result_str.split('\n')])
    aln_df = pd.read_csv(io.StringIO(align_result_str),
        sep='\t',
        header=None, # blastn returns no header in "outfmt 6"
        names= ['qseqid', 'sseqid', 'sstrand',
                'length',   'qlen',    'slen',
                'qstart',   'qend',  'sstart',
                                       'send'],
        true_values=['plus'], # plus strand will be True
        false_values=['minus'] # minus strand will be False
    )

    return aln_df
# end def str2df


def make_outfpath(fq_fpath):
    # Function configures path to output file.
    #
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath: str;
    #
    # Returns path to output file (str).

    return os.path.join(
        os.path.dirname(fq_fpath),
        'cleaned_{}'.format(os.path.basename(fq_fpath).replace('.gz', ''))
    )
# end def make_outfpath


def fastq_records(fq_fpath):
    # Function (generator) that "generates" fastq records from given file.
    #
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath: str;

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
    # Function writes a fastq record to given file
    #
    # :param fq_record: dictionary of following structure:
    #   {'seq_id': <seq_title>, 'seq': <sequence>, 'cmnt': <comment>, 'qual': <quality_string>};
    # :type fq_record: dict<str: str>;
    # :param outfile: file to which a record will be written;
    # :type outfile: _io.TextIOWrapper;

    outfile.write('@{}\n{}\n{}\n{}\n'.format(fq_record['seq_id'],
            fq_record['seq'], fq_record['cmnt'], fq_record['qual'])
    )
# end def write_fastq_record


# Positions in BLAST+ are 1-based
MAX_EDGE_OFFSET = 1


def set_touch_start(row):
    # Function to be used with pandas.DataFrame.apply function.
    # It adds a column indicating if the alignment in a particular row
    #   "touches" 5'-end of subject sequence.
    #
    # :param row: row of dataframe to which this function if applied;
    # :type row: pandas.Series;

    if row['sstrand']:
        row['touch_start'] = row['sstart'] < MAX_EDGE_OFFSET + 2
    else:
        row['touch_start'] = row['sstart'] > row['slen'] - 2
    # end if

    return row
# end def get_touch_start


def set_touch_end(row):
    # Function to be used with pandas.DataFrame.apply function.
    # It adds a column indicating if the alignment in a particular row
    #   "touches" 3'-end of subject sequence.
    #
    # :param row: row of dataframe to which this function if applied;
    # :type row: pandas.Series;

    if row['sstrand']:
        row['touch_end'] = row['send'] > row['slen'] - 1 - MAX_EDGE_OFFSET
    else:
        row['touch_end'] = row['send'] < MAX_EDGE_OFFSET + 2
    # end if

    return row
# end def get_touch_start


def get_aligned_spans(curr_alns):
    # Function analyses obtained alignments and find "spans" -- subsequences
    #   into which input read should be "shredded".
    # These spans should not overlap
    #
    # :param curr_alns: dataframe containing alignmens data obtained from str2df function;
    # :type curr_alns: pandas.DataFrame;

    # If there are no alignments -- just return empty list.
    if curr_alns.empty:
        return []
    # end if

    # Add columns indicating if alignments "touch" 5'- and 3'-ends of subjects
    curr_alns = curr_alns.apply(set_touch_start, axis=1)
    curr_alns = curr_alns.apply(set_touch_end, axis=1)

    # Extract "touching" alignments
    full_span_alns = curr_alns.query('touch_start & touch_end') # from 5' to 3'
    one_side_alns = curr_alns.query('@operator.xor(touch_start, touch_end)') # from 5' of from 3' (with a break somewhere in between)

    # List of result spans
    aligned_spans = list()

    # Get major fragment first
    # Their accessions start with 'A' and will appear firstly in `full_span_alns`
    full_span_alns = full_span_alns.sort_values(by='sseqid', ascending=True)
    # full_span_alns.sort_values(by='sseqid', ascending=True, inplace=True)

    # We will use this array to make sure that one span
    #   is included in result ("shredded") reads multiple times.
    cov_array = np.zeros(curr_alns.iloc[0,:]['qlen'], dtype=np.uint8)

    # Find spans
    # Firsly check full_span_alns, since they are major product and should preceed
    for aln_collenstion in (full_span_alns, one_side_alns):

        for i in range(aln_collenstion.shape[0]):

            aln = aln_collenstion.iloc[i, :] # get current alignment
            buff_cov_array = np.copy(cov_array) # copy current `cov_array` to buffer
            np.add.at(buff_cov_array, range(aln['qstart']-1, aln['qend']), 1) # add 1-s to find overlaps

            # If it is a 2 (number two) somewhere, we have overlap and do not need this span.
            # if no 2 emerger, keep this span and add it to `aligned_spans`
            if not any(buff_cov_array > 1):
                np.add.at(cov_array, range(aln['qstart']-1, aln['qend']), 1) # update cov_array
                aligned_spans.append( (aln['qstart']-1, aln['qend']) )
            # end if
        # end for
    # end for

    return aligned_spans
# end def get_aligned_spans


def clean_and_shred(aln_df, fq_fpath):
    # Function organizes analysis of dataframe of alignment data.
    #
    # :param aln_df: dataframe containing alignment data to be analysed;
    # :type aln_df: pandas.DatFrame;
    # :param fq_fpath: path to input fastq file;
    # :type fq_fpath str;
    #
    # Returns path to output file.

    # Get path to output file.
    outfpath = make_outfpath(fq_fpath)

    # Count reads and configure variables for printing status bar
    nreads = sum(1 for _ in OPEN_FUNCS[int(fq_fpath.endswith('.gz'))](fq_fpath)) // 4
    bar_len = int(os.get_terminal_size().columns * 0.50)
    next_print_num = int(nreads * 0.01)
    inc_num = next_print_num

    print('{} - Start cleaning and shredding...'.format(getwt()))

    sys.stdout.write('{} - [{}] 0/{}'.format(getwt(), ' '*bar_len, nreads))
    sys.stdout.flush()

    with open(outfpath, 'w') as outfile:

        for i, fq_record in enumerate(fastq_records(fq_fpath)):

            # Select rows containing alignments of current read
            # curr_alns = aln_df[aln_df['qseqid'] == fq_record['seq_id']]
            curr_alns = aln_df.query('qseqid == @fq_record["seq_id"]')

            # Analyse alignments and get spans
            aligned_spans = get_aligned_spans(curr_alns)

            # Write spans
            for aln_span in aligned_spans:

                # Configure variables for output fastq record
                qstart, qend = aln_span[0], aln_span[1]
                curr_seq_id = '{}_{}-{}'.format(fq_record['seq_id'], qstart+1, qend) # write 1-based coordinates here
                curr_seq = fq_record['seq'][qstart : qend]
                curr_qual = fq_record['qual'][qstart : qend]

                write_fastq_record({'seq_id':  curr_seq_id,
                                    'seq'   :  curr_seq,
                                    'cmnt'  :  fq_record['cmnt'],
                                    'qual'  :  curr_qual},
                                   outfile)
            # end for

            # Update status bar
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


def rm_query_file(query_fpath):
    # Function removes temporary query file.
    #
    # :param query_fpath: path to query file;
    # :type query_fpath: str;

    try:
        os.unlink(query_fpath)
    except OSError as oserr:
        print('Warning: Cannot remove temporary file `{}`'.format(query_fpath))
        print( str(oserr) )
    # end try
# end rm_query_file


def main():
    # Handle command line arguments
    fq_fpaths, db_fpath, n_thr = handle_cl_args()

    print('{} - Start.'.format(getwt()))

    for fq_fpath in fq_fpaths:
        print('{} - Processing file `{}`'.format(getwt(), fq_fpath))

        # Convert input fastq file to fastq format in order to pass the latter to blastn
        query_fpath = fastq2fasta(fq_fpath)

        # Align and obtain dataframe containing data about alignments
        aln_df = str2df(disco_align(query_fpath, db_fpath, n_thr))

        # Analyse alignments, configure spans and write them to output file in fastq format
        outfpath = clean_and_shred(aln_df, fq_fpath)

        # Remove temporary query file
        rm_query_file(query_fpath)

        print('{} - File `{}` is processed.'.format(getwt(), fq_fpath))
        print('Output file: `{}`'.format(outfpath))
        print('-------')
    # end for
# end def main


if __name__ == '__main__':
    main()
# end if