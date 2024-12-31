#!/usr/bin/env python3

# The script prints which portion of a read is actually mapped,
#   and which portions are clipped from each side.
# It might be useful for manual inspection of structural variants.
# Example:
#  Input (a CIGAR string):
#    10H136M1I341M1D8M2D725M23912H
#  Output:
#    10-H|          1,211|       23,912-H

# Usage:
# Mode 1: pass CIGAR strings as arguments
#   $ ./calc_cigar_qlen.py 32S67M8S
#   $ ./calc_cigar_qlen.py 32S67M8S 12S67M1D34M18S
# Mode 2: read CIGAR lines from stdin
#   $ samtools view my.bam | grep 'READ_ID' | cut -f6 | ./calc_cigar_qlen.py

# __version__ = '1.0.a'


import re
import sys


CIGAR_PATTERN = re.compile(
    r'^[0-9MIDNSHPX=]+$'
)


def main():
    if len(sys.argv) > 1:
        # Mode 1: parse CIGAR from arguments passed
        for cigar_str in sys.argv[1:]:
            validate_cigar(cigar_str)
            print_cigar_qlen(cigar_str)
        # end for
    else:
        # Mode 2: read CIGAR lines from stdin
        eof = False
        while not eof:
            line = sys.stdin.readline().strip()
            if line != '':
                validate_cigar(line)
                print_cigar_qlen(line)
            else:
                eof = True
            # end if
        # end while
    # end if
# end if


def validate_cigar(cigar_str):
    rematch = re.match(CIGAR_PATTERN, cigar_str)
    if rematch is None:
        sys.stderr.write(
            'Error: invalid CIGAR: `{}`'.format(cigar_str)
        )
        sys.exit(1)
    # end if
# end def


def print_cigar_qlen(cigar_str):
    q_cigar_operations = re.findall(r'[0-9]+M', cigar_str) \
                       + re.findall(r'[0-9]+X', cigar_str) \
                       + re.findall(r'[0-9]+=', cigar_str) \
                       + re.findall(r'[0-9]+I', cigar_str)

    lengths = map(
        lambda x: int(x[:-1]),
        q_cigar_operations
    )
    q_len = '{:,}'.format(sum(lengths))

    start_clip_rematch = re.search(r'^([0-9]+[HS])', cigar_str)
    if not start_clip_rematch is None:
        clip_str = start_clip_rematch.group(1)
        start_clip = '{:,}-{}'.format(int(clip_str[:-1]), clip_str[-1])
    else:
       start_clip = '.'
    # end if
    end_clip_rematch = re.search(r'([0-9]+[HS])$', cigar_str)
    if not end_clip_rematch is None:
        clip_str = end_clip_rematch.group(1)
        end_clip = '{:,}-{}'.format(int(clip_str[:-1]), clip_str[-1])
    else:
        end_clip = '.'
    # end if

    sys.stdout.write(
        '{:>15}|{:>15}|{:>15}\n'.format(start_clip, q_len, end_clip)
    )
# end def


if __name__ == '__main__':
    main()
# end if
