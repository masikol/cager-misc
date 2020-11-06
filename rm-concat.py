#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.1.a"
# Year, month, day
__last_update_date__ = "2020-11-05"

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
import multiprocessing as mp


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


rc_dict = {
    'A': 'T', 'T': 'A',
    'G': 'C', 'C': 'G',
    'N': 'N'
}

def single_nucl_rc(nucl):
    return rc_dict[nucl]
# end def single_nucl_rc

def rc(seq):
    return "".join(map(single_nucl_rc, seq[::-1]))
# end def rc


# -- Smith-Waterman --

from array import array
from functools import partial

class SWAlignLensNotEq(Exception):
    """
    An exception meant to be raised when length of the first aligned 
    sequence the second's one are not equal after aligning.
    Is subclass of Exception.
    """
    pass
# end class SWAlignLensNotEq


class SWInvalidMove(Exception):
    """
    An exception meant to be raised when Smith-Waterman algorithm can't perform the next move.
    Is subclass of Exception.
    """
    pass
# end class SWInvalidMove


class AlignResult:
    """
    Class AlignResult is dedicated to perform result of Smith-Waterman aligning.

    :field q_align: aligned first (aka query) sequence (e.i. with gaps represented as '-');
    :type q_align: str;
    :field s_align: aligned seconf (aka subject) sequence (e.i. with gaps represented as '-');
    :type s_align: str;
    :field q_start: 1-based number of query sequence, in which alignment starts;
    :type q_start: int;
    :field q_end: 1-based number of query sequence, in which alignment ends;
    :type q_end: int;
    :field s_start: 1-based number of subject sequence, in which alignment starts;
    :type s_start: int;
    :field s_end: 1-based number of subject sequence, in which alignment ends;
    :type s_end: int;
    """

    def __init__(self, q_align, s_align, q_start, q_end, s_start, s_end, qlen,
        score, matrix):
        """
        :param q_align: aligned first (aka query) sequence (e.i. with gaps represented as '-');
        :type q_align: str;
        :param s_align: aligned seconf (aka subject) sequence (e.i. with gaps represented as '-');
        :type s_align: str;
        :param q_start: 1-based number of query sequence, in which alignment starts;
        :type q_start: int;
        :param q_end: 1-based number of query sequence, in which alignment ends;
        :type q_end: int;
        :param s_start: 1-based number of subject sequence, in which alignment starts;
        :type s_start: int;
        :param s_end: 1-based number of subject sequence, in which alignment ends;
        :type s_end: int;
        """

        self.q_align = q_align
        self.s_align = s_align

        self.q_start = q_start
        self.q_end = q_end
        self.s_start = s_start
        self.s_end = s_end

        self.qlen = qlen

        self.score= score
        self.matrix = matrix
    # end def __init__

    def __repr__(self):
        return "\nQS={}; QE={}; SS={}; SE={}".format( self.q_start, self.q_end,
                                                            self.s_start, self.s_end )
    # end def __repr__
# end class AlignResult


def SW_get_new_score(up, left, diag,
                  matched, gap_penalty, match, mismatch):
    """
    Function is dedicated to fill an element of alignment matrix during
    forward step of Smith-Waterman algorithm.

    :param up: value of the upper element of matrix;
    :type up: int;
    :param left: value of the left element of matrix;
    :type left: int;
    :param diag:value of the diagonal element of matrix;
    :type diagonal: int;
    :param matched: logic value: is True if nucleotides are equal and False otherwise;
    :type matched: bool;
    :param gap_penalty: penalty for gap;
    :type gap_penalty: int;
    :param match: reward for matching;
    :type match: int;
    :param mismatch: penalty for mismatch;
    :type mismatch: int;
    """

    add_to_diag = match if matched else mismatch
    return max(0, diag+add_to_diag, left-gap_penalty, up-gap_penalty)
# end def SW_get_new_score


# Constants for traceback.
END, DIAG, UP, LEFT = range(4)


def SW_next_move(SWsm, i, j, match, mismatch, gap_penalty, matched):
    """
    Function is dedicated to perform a move during traceback of Smith-Waterman algorithm.

    :param SWsm: Smith-Waterman scoring matrix;
    :type SWsm: list<list<int>>;   meant to be changed to list<array.array<unsigned_int>>
    :param i: row index of scoring matrix;
    :type i: int;
    :param j: column index of scoring matrix;
    :type j: int;
    :param gap_penalty: penalty for gap;
    :type gap_penalty: int;
    :param matched: logic value: is True if nucleotides are equal and False otherwise;
    :type matched: bool;

    Returns 1 if move is diagonal, 2 if move is upper, 3 if move is left.

    Raises SWInvalidMove when Smith-Waterman algorithm can't perform the next move.
    """
    diag = SWsm[i - 1][j - 1]
    up   = SWsm[i - 1][  j  ]
    left = SWsm[  i  ][j - 1]

    if (diag + (mismatch, match)[matched]) == SWsm[i][j]:
        return DIAG if diag!=0 else END
    
    elif (up - gap_penalty) == SWsm[i][j]:
        return UP if up!=0 else END
    
    elif (left - gap_penalty) == SWsm[i][j]:
        return LEFT if up!=0 else END
    
    else:
        # Execution should not reach here.
        raise SWInvalidMove("""Smith-Waterman algorithm: invalid move during traceback:
            diag={}; up={}; left={}; matched={}; SWsm[i][j]={}""".format(diag, up, left, matched, SWsm[i][j]))
    # end if
# end def SW_next_move


def SW_align(jseq, iseq, gap_penalty=10, match=1, mismatch=-1):
    """
    Function performs Smith-Waterman aligning algorithm.
    Reward and all penalties are int for the sake of acceleration of this slow snake.

    :param jseq: 'query' sequence, which is loated to the top of scoring matrix, across the j index varying;
    :type jseq: str;
    :param iseq: 'subject' sequence, which is loated to the left of scoring matrix, across the i index varying;
    :type iseq: str;
    :param gap_penalty: penalty for gap;
    :type gap_penalty: int;
    :param match: reward for matching;
    :type match: int;
    :param mismatch: penalty for mismatch;
    :type mismatch: int;

    Returns instance of AlignResult performing results of alignment.
    """

    SWsm = [array('I', [0 for j in range(len(jseq)+1)]) for i in range(len(iseq)+1)]
    # SWsm = [[0 for j in range(len(jseq)+1)] for i in range(len(iseq)+1)]
    score = partial(SW_get_new_score, gap_penalty=gap_penalty, match=match, mismatch=mismatch)
    max_val, max_pos = 0, None

    for i in range(1, len(iseq)+1):
        for j in range(1, len(jseq)+1):
            SWsm[i][j] = score(up=SWsm[i-1][j], left=SWsm[i][j-1], diag=SWsm[i-1][j-1], matched=iseq[i-1] == jseq[j-1])

            if SWsm[i][j] > max_val:
                max_val = SWsm[i][j]
                max_pos = (i, j)
            # end if
        # end for
    # end for

    # It happens (no alignment)
    if max_pos is None:
        return None
    # end if

    aligned_jseq = ""
    aligned_iseq = ""
    i, j = max_pos
    move = SW_next_move(SWsm, i, j, match, mismatch,
                        gap_penalty, iseq[i-1] == jseq[j-1])

    while move != END:

        if move == DIAG:
            aligned_jseq += jseq[j-1]
            aligned_iseq += iseq[i-1]
            j -= 1
            i -= 1

        elif move == UP:
            aligned_jseq += '-'
            aligned_iseq += iseq[i-1]
            i -= 1

        elif move == LEFT:
            aligned_jseq += jseq[j-1]
            aligned_iseq += '-'
            j -= 1
        # end if

        move = SW_next_move(SWsm, i, j, match, mismatch,
                            gap_penalty, iseq[i-1] == jseq[j-1])
    # end while

    if jseq[j-1] == iseq[i-1]:
        aligned_jseq += jseq[j-1]
        aligned_iseq += iseq[i-1]
    # end if

    return AlignResult(aligned_jseq[::-1], aligned_iseq[::-1], j, max_pos[1], i-1, max_pos[0], len(jseq),
        max_val, SWsm)
# end def SW_align

# --------------------


def handle_cl_args():

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'p:t:', ['primers=', 'threads='])
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
    n_thr = 1
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

        elif opt in ("-t", "--threads"):
            try:
                n_thr = int(arg)
                if n_thr < 1:
                    raise ValueError
                # end if
            except ValueError:
                print_error("number of threads must be positive integer number!")
                print(" And here is your value: '{}'".format(arg))
                sys.exit(1)
            # end try

            if n_thr > len(os.sched_getaffinity(0)):
                print("""\nWarning! You have specified {} threads to use
      although {} are available.""".format(n_thr, len(os.sched_getaffinity(0))))
                error = True
                while error:
                    reply = input("""\nPress ENTER to switch to {} threads,
      or enter 'c' to continue with {} threads,
      or enter 'q' to exit:>>""".format(len(os.sched_getaffinity(0)), n_thr))
                    if reply in ("", 'c', 'q'):
                        error = False
                        if reply == "":
                            n_thr = len(os.sched_getaffinity(0))
                            print("\nNumber of threads switched to {}\n".format(n_thr))
                        elif reply == 'c':
                            pass
                        elif reply == 'q':
                            sys.exit(0)
                        # end if
                    else:
                        print("\nInvalid reply!\n")
                    # end if
                # end while
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
        primers_seqs = list(map(str.upper,                # to upper case
            map(str.strip,                                 # srtip
                filter(is_seq, primers_file.readlines())   # extract only sequences
               )
            )
        )
    # end with

    for i in range(1, len(primers_seqs), 2):
        primers_seqs[i] = rc(primers_seqs[i])
    # end for
    flow_range = [i-1 if i%2 == 0 else i+1 for i in range(1, len(primers_seqs)-1)]
    flow_range.insert(0, 0)
    flow_range.append(len(primers_seqs)-1)

    primers_flow = [primers_seqs[i] for i in flow_range]

    return fq_fpaths, primers_flow, n_thr
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


def write_fastq_record(fq_record, outfpath):
    with open(outfpath, 'a') as outfile:
        outfile.write('{}\n{}\n{}\n{}\n'.format(fq_record['seq_id'],
                fq_record['seq'], fq_record['cmnt'], fq_record['qual'])
        )
    # end with
# end def write_fastq_record


def check_repeats(read_seq, primers_seqs):
    is_concat = False
    for seq in primers_seqs:
        if read_seq.count(seq) > 1:
            is_concat = True
            break
        # end if
    # end for
    return is_concat
# end def check_repeats


def check_order(fq_record, primers_seqs):

    is_concat = True
    read_seq = fq_record['seq']

    first_match = False
    start_fragment_i = None
    start_fragment_pos = None
    max_perm_mismatches = 2

    for i, prim_seq in enumerate(primers_seqs):

        min_perm_score = len(prim_seq) - max_perm_mismatches
        aligmnent = SW_align(prim_seq, read_seq)

        if aligmnent.score >= min_perm_score:
            first_match = True
            start_fragment_i = i
            start_fragment_pos = aligmnent.s_start
            break
        # end if
    # end for

    if first_match and start_fragment_i < (len(primers_seqs)-1):

        prim_seq = primers_seqs[start_fragment_i+1]
        min_perm_score = len(prim_seq) - max_perm_mismatches
        aligmnent = SW_align(prim_seq, read_seq)

        if (aligmnent.score >= min_perm_score) & (aligmnent.s_start > start_fragment_pos):

            end_fragment_i = start_fragment_i+1
            end_fragment_pos = aligmnent.s_end

            is_concat = False

            for i in (2, 3):
                try:
                    prim_seq = primers_seqs[start_fragment_i+i]
                except IndexError:
                    break
                # end try
                min_perm_score = len(prim_seq) - max_perm_mismatches
                aligmnent = SW_align(prim_seq, read_seq)

                if aligmnent.score >= min_perm_score & aligmnent.s_start > start_fragment_pos:
                    end_fragment_i += 1
                    end_fragment_pos = aligmnent.s_end
                # end if
            # end for
        # end if
    # end if

    if is_concat:
        return is_concat, None
    else:
        fq_record['seq'] = read_seq[start_fragment_pos : end_fragment_pos]
        fq_record['qual'] = fq_record['qual'][start_fragment_pos : end_fragment_pos]
        return is_concat, fq_record
    # end if
# end def check_order


def concat_clean(data, fq_fpath, primers_seqs, nreads, ok_out_fpath, trash_out_fpath):

    bar_len = int(os.get_terminal_size().columns * 0.50) - 1 # minus one for >

    rc_primers_seqs = tuple(map(rc, primers_seqs))[::-1]

    for fq_record in data:

        read_seq = fq_record['seq'].upper()
        is_concat = False
        
        forw_repeats = check_repeats(read_seq, primers_seqs)
        rc_repeats = check_repeats(read_seq, rc_primers_seqs)
        is_concat = forw_repeats | rc_repeats

        if not is_concat:
            is_concat, trimmed_fq_record = check_order(fq_record, primers_seqs)
            if is_concat:
                is_concat, trimmed_fq_record = check_order(fq_record, rc_primers_seqs)
            # end if
        # end if

        with write_lock:
            if is_concat:
                write_fastq_record(fq_record, trash_out_fpath)
                trash_count.value += 1
            else:
                write_fastq_record(trimmed_fq_record, ok_out_fpath)
                ok_count.value += 1
            # end if
        # end with

        with counter_lock:
            counter.value += 1
            done_i = counter.value
        # end with

        done_ratio = done_i / nreads
        with print_lock:
            sys.stdout.write('\r[{}>{}] {}/{}'.format('='*int(bar_len*done_ratio),
                ' '*int(bar_len*(1-done_ratio)), done_i, nreads) )
            sys.stdout.flush()
        # end with
    # end for

    return ok_count, trash_count
# end def concat_clean


def proc_init(write_lock_buff, print_lock_buff, counter_buff, counter_lock_buff, ok_count_buff, trash_count_buff):

    global write_lock
    write_lock = write_lock_buff

    global print_lock
    print_lock = print_lock_buff

    global counter
    counter = counter_buff

    global counter_lock
    counter_lock = counter_lock_buff

    global ok_count
    ok_count = ok_count_buff

    global trash_count
    trash_count = trash_count_buff
# end def proc_init


def read_fastq_record(fq_file, fmt_func):

    fastq_rec = {                    #read all 4 lines of fastq-record
        "seq_id": fmt_func(fq_file.readline()),
        "seq": fmt_func(fq_file.readline()).upper(), # searching for cross-talks is case-dependent
        "cmnt": fmt_func(fq_file.readline()),
        "qual": fmt_func(fq_file.readline())
    }

    # If all lines from files are read
    if fastq_rec["seq_id"] == "":
        return None
    # end if
    return fastq_rec
# end def read_fastq_record


def fastq_read_packets(fq_fpath, nreads, n_thr):

    how_to_open = OPEN_FUNCS[ fq_fpath.endswith('.gz') ]
    fmt_func = FORMATTING_FUNCS[ fq_fpath.endswith('.gz') ]

    fq_file = how_to_open(fq_fpath)

    # Compute packet size (one packet -- one thread).
    pack_size = nreads // n_thr
    if nreads % n_thr != 0:
        pack_size += 1
    # end if

    for i in range(n_thr):
        packet = list()
        for j in range(pack_size):
            tmp = read_fastq_record(fq_file, fmt_func)
            # if the end of the line is reached
            if tmp is None:
                # Yield partial packet and return 'None' on next call
                yield packet
                return
            else:
                packet.append(tmp)
            # end if
        # end for
        yield packet # yield full packet
    # end for
# end def fastq_read_packets


def main():
    fq_fpaths, primers_seqs, n_thr = handle_cl_args()

    ok_count, trash_count = 0, 0
    print('{} - Start.'.format(getwt()))

    for fq_fpath in fq_fpaths:
        print('{} - Processing file `{}`.'.format(getwt(), fq_fpath))

        bar_len = int(os.get_terminal_size().columns * 0.50) - 1 # minus one for >
        nreads = sum(1 for _ in OPEN_FUNCS[int(fq_fpath.endswith('.gz'))](fq_fpath)) // 4

        sys.stdout.write('[{}] 0/{}'.format(' '*bar_len, nreads))
        sys.stdout.flush()

        ok_out_fpath, trash_out_fpath = make_outfpaths(fq_fpath)
        for f in (ok_out_fpath, trash_out_fpath):
            with open(f, 'w') as file:
                pass
            # end with
        # end for

        ok_count_tmp = mp.Value('i', 0)
        trash_count_tmp = mp.Value('i', 0)

        pool = mp.Pool(n_thr, initializer=proc_init,
            initargs=(mp.Lock(), mp.Lock(), mp.Value('i', 0), mp.Lock(),
                ok_count_tmp, trash_count_tmp))
        pool.starmap(concat_clean,
            [(data, fq_fpath, primers_seqs, nreads, ok_out_fpath, trash_out_fpath)
            for data in fastq_read_packets(fq_fpath, nreads, n_thr)])

        # Reaping zombies
        pool.close()
        pool.join()

        ok_count += ok_count_tmp.value
        trash_count += trash_count_tmp.value

        sys.stdout.write('\r[{}] {}/{}\n'.format('='*bar_len, nreads, nreads))
        sys.stdout.flush()

        print('{} - File `{}` is procesed.'.format(getwt(), fq_fpath))
        print('Result files:')
        for f in (ok_out_fpath, trash_out_fpath):
            print(' `{}`'.format(f))
        # end for
        print('-'*10)
    # end for

    print('Results: {} reads keeped, {} reads discarded.'.format('{:,}'.format(ok_count).replace(',', ' '),
                                                                 '{:,}'.format(trash_count).replace(',', ' '))
    )
# end def main


if __name__ == '__main__':
    main()
# end if