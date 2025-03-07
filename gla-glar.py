#!/usr/bin/env python3

import os
import re
import sys
import argparse
import statistics as sts
from typing import Sequence


__version__ = '1.0.a'
__author__  = 'Maksim Sikolenko'
__last_update_date__ = '2025-03-07'

sys.stderr.write(
    '# INFO: gla-glar.py: GLAde GLARer version {}, {}, by {}\n'.format(
        __version__, __last_update_date__, __author__
    )
)
sys.stderr.flush()


# >>> Parse arguments >>>

parser = argparse.ArgumentParser()

parser.add_argument(
    'input_fpath',
    help='Input file in Genbank format'
)

parser.add_argument(
    '-g',
    '--gap-length',
    help='Minimum length of an integenic region to report',
    type=int,
    default=700
)

parser.add_argument(
    '-m',
    '--mean-frac',
    help='''Intergenic regions longer than the specified value*mean will be reported, for each sequence.
For example `-m 2` will report regions with length 2*mean.''',
    type=float
)

parser.add_argument(
    '-M',
    '--median-frac',
    help='''Intergenic regions longer than the specified value*median will be reported, for each sequence.
For example `-M 2` will report regions with length 2*median.''',
    type=float
)


args = parser.parse_args()


input_fpath = os.path.abspath(args.input_fpath)
if not os.path.exists(input_fpath):
    sys.stderr.write('Error! Input file does not exist: `{}`!\n'.format(input_fpath))
    sys.exit(1)
# end if

min_gap_len = None
min_gap_len_times_mean = None
min_gap_len_times_median = None
if not args.gap_length is None:
    min_gap_len = int(args.gap_length)
# end if
if not args.mean_frac is None:
    min_gap_len_times_mean = float(args.mean_frac)
# end if
if not args.median_frac is None:
    min_gap_len_times_median = float(args.median_frac)
# end if

for x in min_gap_len, min_gap_len_times_mean, min_gap_len_times_median:
    if not x is None and x < 1:
        sys.stderr.write('Error: values < 1 are not allowed: `{}`\n'.format(x))
        sys.exit(1)
    # end if
# end for
del x

mean_mode = not min_gap_len_times_mean is None
median_mode = not min_gap_len_times_median is None
if mean_mode and median_mode:
    sys.stderr.write('Error: please specify either `-m` or `-M`, but not both.`\n')
    sys.exit(1)
# end if

mean_median_mode = mean_mode or median_mode


# <<< Parse arguments <<<


# >>> Functions >>>

def parse_feature_coordinates(record_str : str) -> Sequence[Sequence[int]]:
    feature_intervals = set()

    gene_pattern = re.compile(
        r'^[ ]*gene[ ]*(complement\()?([0-9]+)\.\.([0-9]+)(\))?$'
    )
    cds_pattern = re.compile(
        r'^[ ]*CDS[ ]*(complement\()?([0-9]+)\.\.([0-9]+)(\))?$'
    )

    for line in record_str.splitlines():
        line = line.strip()
        gene_match = re.match(gene_pattern, line)
        cds_match  = re.match(cds_pattern,  line)
        for match_obj in (gene_match, cds_match):
            if not match_obj is None:
                feature_intervals.add(
                    parse_interval(match_obj)
                )
            # end if
        # end for
    # end for

    return tuple(
        sorted(
            feature_intervals,
            key=lambda x: x[0]
        )
    )
# end def

def parse_interval(match_obj : re.Match) -> Sequence[int]:
    left_coord  = int(match_obj.group(2))
    right_coord = int(match_obj.group(3))
    return left_coord, right_coord
# end def


def make_interfeature_intervals(feature_intervals : Sequence[Sequence[int]]) -> Sequence[Sequence[int]]:
    num_interfeature_intervals = len(feature_intervals) - 1

    interfeature_intervals = [
        None for i in range(num_interfeature_intervals)
    ]

    for i in range(num_interfeature_intervals):
        interfeature_intervals[i] = \
            (
                feature_intervals[ i ][1] - 1,
                feature_intervals[i+1][0] + 1,
            )
    # end for

    return interfeature_intervals
# end def


def make_interfeature_lengths(interfeature_intervals : Sequence[Sequence[int]]) -> Sequence[int]:
    return tuple(
        map(
            lambda interval: interval[1] - interval[0] + 1,
            interfeature_intervals
        )
    )
# end def


def print_long_intervals_simple(interfeature_intervals : Sequence[Sequence[int]],
                                min_gap_len : int):

    threshold = min_gap_len - 1 # > is expected to be faster than >=
    sys.stderr.write('# INFO: Reporting intergenic regions longer than {} bp.\n'.format(threshold))
    sys.stderr.flush()

    for interval in interfeature_intervals:
        interval_length = interval[1] - interval[0] + 1
        if interval_length > threshold:
            print('{}-{}'.format(interval[0], interval[1]))
        # end if
    # end for
# end def


def print_long_intervals_mM(interfeature_intervals : Sequence[Sequence[int]],
                            interfeature_lengths : Sequence[int],
                            min_gap_len_times_mean : float,
                            min_gap_len_times_median : float):
    if not min_gap_len_times_mean is None:
        aggregate_fun = sts.mean
        frac = min_gap_len_times_mean
    elif not min_gap_len_times_median is None:
        aggregate_fun = sts.median
        frac = min_gap_len_times_median
    # end if

    mean_or_median = round(aggregate_fun(interfeature_lengths))

    if not min_gap_len_times_mean is None:
        sys.stderr.write('# INFO: Mean length of an intergenic region: {} bp.\n'.format(mean_or_median))
    elif not min_gap_len_times_median is None:
        sys.stderr.write('# INFO: Median length of an intergenic region: {} bp.\n'.format(mean_or_median))
    # end if

    threshold = round(frac * mean_or_median) - 1 # > is expected to be faster than >=
    sys.stderr.write('# INFO: Reporting intergenic regions longer than {} bp.\n'.format(threshold))
    sys.stderr.flush()

    for i, interval_len in enumerate(interfeature_lengths[:-1]):
        if interval_len > threshold:
            print(
                '{}-{}'.format(
                    interfeature_intervals[i][1],
                    interfeature_intervals[i+1][0],
                )
            )
        # end if
    # end for
# end def

# <<< Functions <<<


# >>> Proceed >>>

with open(input_fpath, 'rt') as input_handle:
    input_records = tuple(
        map(
            str.strip,
            input_handle.read().split('\n//\n')
        )
    )
# end with

input_records = tuple(
    filter(
        lambda x: x != '',
        input_records
    )
)


for record_str in input_records:
    first_line = record_str.partition('\n')[0]
    sys.stdout.write('{}\n'.format(first_line))

    feature_intervals = parse_feature_coordinates(record_str)

    interfeature_intervals = make_interfeature_intervals(
        feature_intervals
    )
    del feature_intervals

    if not mean_median_mode:
        print_long_intervals_simple(interfeature_intervals, min_gap_len)
    else:
        interfeature_lengths = make_interfeature_lengths(interfeature_intervals)
        print_long_intervals_mM(
            interfeature_intervals,
            interfeature_lengths,
            min_gap_len_times_mean,
            min_gap_len_times_median
        )
    # end if
# end for

sys.exit(0)
