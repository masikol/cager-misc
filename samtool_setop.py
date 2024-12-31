#!/usr/bin/env python3

# The script performs set operations (intersection, unuin, difference) on read IDs that are mapped to different reference positions.

# Usage:
# python3 samtools_setop.py [-f FLAG] [-F FLAG] <BAM> <RNAME> <POS1> <SET_OPERATION> <POS2>
# allowed SET_OPERATION values:
#   'I': Intersection
#   'D': Difference
#   'U': Union

# Example:
# python3 samtools_setop.py \
#     v1b_nan_F4.sorted.bam \
#     Lpb_B-492_v1b \
#     463585 I 469240

# __version__ = '1.0.a'


import os
import sys
import argparse
import subprocess as sp


def main():
    args = parse_args()
    validate_args(args)
    fF_flags_str = make_fF_flags_str(args.f, args.F)
    read_ids_1 = get_read_ids_for_pos(args.bam, args.rname, args.pos_1, fF_flags_str)
    read_ids_2 = get_read_ids_for_pos(args.bam, args.rname, args.pos_2, fF_flags_str)
    result_read_ids = perform_operation(read_ids_1, read_ids_2, args.operation)
    output_result(result_read_ids)
# end def


def parse_args():
    parser = argparse.ArgumentParser(
        prog='samtools_setop.py',
        description='Set operations on read ids for regions of BAM/SAM files.',
        epilog='Vive valeque!'
    )
    parser.add_argument(
        'bam',
        help='BAM/SAM file of interest.'
    )
    parser.add_argument(
        'rname',
        help='reference name, the same as RNAME column in BAM/SAM files.'
    )
    parser.add_argument(
        'pos_1',
        help='a position within RNAME, left argument of the operation.',
        type=int
    )
    parser.add_argument(
        'operation',
        help='''set operation to perform. One of the following:
                `I`: intersection;
                `U`: union;
                `D`: difference.'''
    )
    parser.add_argument(
        'pos_2',
        help='a position within RNAME, right argument of the operation.',
        type=int
    )
    parser.add_argument(
        '-f',
        help='Consider only reads having all of the FLAGs present.',
        type=int,
        required=False,
        action='append',
        default=list()
    )
    parser.add_argument(
        '-F',
        help='Consider only reads having none of the FLAGs present.',
        type=int,
        required=False,
        action='append',
        default=list()
    )
    args = parser.parse_args()

    args.bam = os.path.abspath(args.bam)
    args.operation = args.operation.upper()

    return args
# end def


def validate_args(args):

    if not os.path.isfile(args.bam):
        sys.stderr.write(
            'Error: input file does not exist: `{}`\n'.format(args.bam)
        )
        sys.exit(1)
    # end if

    allowed_operations = frozenset({
        'I', 'U', 'D'
    })
    if not args.operation in allowed_operations:
        sys.stderr.write(
            'Error: invalid set operation: `{}`\n'.format(args.operation)
        )
        sys.stderr.write(
            'Allowed operations: {}\n'.format(', '.join(allowed_operations))
        )
        sys.exit(1)
    # end if

    positive_int_zip = zip(
        (args.pos_1, args.pos_2),
        ('pos_1', 'pos_2'),
    )
    for arg, argname in positive_int_zip:
        validate_positive_integer(arg, argname)
    # end for
    for f_value in args.f:
        validate_positive_integer(f_value, 'f')
    # end for
    for F_value in args.F:
        validate_positive_integer(F_value, 'F')
    # end for
# end def

def validate_positive_integer(arg, argname):
    if arg < 1:
        sys.stderr.write(
            'Error: {} must be an integer > 0. Got `{}`\n'.format(
                argname,
                arg
            )
        )
        sys.exit(1)
    # end if
# end def


def make_fF_flags_str(args_f, args_F):
    f_flags_str = ' '.join(
        map(
            lambda flag: '-f {}'.format(flag),
            args_f
        )
    )
    F_flags_str = ' '.join(
        map(
            lambda flag: '-F {}'.format(flag),
            args_F
        )
    )
    return ' '.join((f_flags_str, F_flags_str))
# end def


def get_read_ids_for_pos(bam_fpath, rname, pos, fF_flags_str):

    cmd = ' '.join(
        [
            'samtools view',
            fF_flags_str,
            bam_fpath,
            '{}:{}-{}'.format(rname, pos, pos),
        ]
    )

    pipe = sp.Popen(
        cmd,
        shell=True,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        encoding='utf-8'
    )
    outs, errs = pipe.communicate()

    if pipe.returncode != 0:
        sys.stderr.write('Error!\n')
        sys.stderr.write(errs)
        sys.stderr.write('\n')
        sys.exit(1)
    # end if

    lines = map(
        lambda l: l.split('\t')[0].strip(),
        outs.splitlines()
    )
    return frozenset(lines)
# end def


def perform_operation(read_ids_1, read_ids_2, operation):
    if operation == 'I':
        return read_ids_1 & read_ids_2
    elif operation == 'U':
        return read_ids_1 | read_ids_2
    elif operation == 'D':
        return read_ids_1 - read_ids_2
    # end if
# end def


def output_result(result_read_ids):
    for s in sorted(result_read_ids):
        print(s)
    # end for
# end def


if __name__ == '__main__':
    main()
# end if
