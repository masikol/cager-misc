#!/usr/bin/env python3
#
# Script moves 'SPAdes-like' *.dna files with coverage less than specified one
#   from 'contigs/' directory to directory 'cov_below_x/'.
# 'SPAdes-like' means that name of file is of following format:
#   NODE_1_length_61704_cov_114.517.dna
# If directory `./contigs/` contains no 'SPAdes-like' files,
#   seqator.py takes them from directory `./contigs/DNA Files`

__version__ = '3.0.a'
__last_update_date__ = '2023-03-01'
__author__ = 'Maxim Sikolenko'
__author_email__ = 'sikolenko@bio.bsu.by'


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
# end if


def platf_depend_exit(exit_code=0):
    """
    Function asks to press ENTER press on Windows
        and exits after that.
    :type exit_code: int;
    """
    if sys.platform.startswith('win'):
        input('Press ENTER to exit:')
    # end if
    sys.exit(exit_code)
# end def

def print_help():
    print('== seqator ==')
    print('Version {}. {} edition.'.format(__version__, __last_update_date__))
    print('By {}, <{}>\n'.format(__author__, __author_email__))

    print('= Description =')
    print('This script performs binning. It bins sequences (in fasta format)')
    print('  and sequence-containing files which have SPAdes-like headers and file names, respectively.')
    print('"SPAdes-like" means the following format:')
    print('  NODE_1_length_61704_cov_114.517')
    print('Seqator can bin according to sequence length and coverage (`length` and `cov` in SPAdes-like line).\n')

    print('= Options =')
    print('-i / --input')
    print('  Input directory or fasta file, depending on seqator_mode (-m).')
    print('  And input fasta file may be gzipped.')
    print('  Mandatory.')
    print('-x / --target-file-extention')
    print('  This option is applicable only if seqator_mode (-m) is `dir`.')
    print('  `-x` is the extention of files to be checked, without the preceding dot.')
    print('  E.g. if you want to bin .fasta files, then specify `-x fasta`.')
    print('  Optional. Default: `dna`.')
    print('# Output')
    print('-o / --output')
    print('  Output directory or fasta file, depending on `-m`.')
    print('  If the file name ends with `.gz`, the output file will be gzipped.')
    print('  Optional. By default, the script will create an output directory in the wokring directory.')
    print('# Seqator mode')
    print('-m / --seqator-mode')
    print('  There are two modes: `dir` and `fasta_file`.')
    print('  If the mode is `dir`, the script will move files which pass the filter')
    print('    from the input directory to the output directory.')
    print('  If the mode is `fasta_file`, the script will copy sequences which pass the filter')
    print('    from the input file to the output file.')
    print('  Also, `-m` may be `auto`. In `auto` mode, the mode will be')
    print('    `dir` if `-i` is a directory and `fasta_file` if `-i` is a regular file.')
    print('  Optional. Default: `auto`.')
    print('# Filter')
    print('-p / --filter-parameter')
    print('  There are two sequence parameters to filter by: `len` and `cov`:')
    print('    length and coverage, respectively.')
    print('  Optional. Default: `cov`.')
    print('-f / --filter-mode')
    print('  There are basically six ways to compare numbers:')
    print('    `lt` (Less Than),    `le` (Less or Equal),')
    print('    `gt` (Greater Than), `ge` (Greater or Equal),')
    print('    `eq` (EQual to),     `ne` (Not Equal to).')
    print('  E.g. if you specify `-p cov -f lt -t 12.5 -m fasta_file`, then')
    print('    the script will copy all sequences having coverage less then 12.5 to the output file.')
    print('  Optional. Default: `lt`.')
    print('-t / --threshold')
    print('  Threshold to use for filtering by `-p` parameter.')
    print('  See the example for `-f` option -- there you\'ll see how `-t` option works.')
    print('  Mandatory.')
    print('# Help and version')
    print('-h / --help')
    print('  Print help message and exit.')
    print('-v / --version')
    print('  Print version and exit.\n')

    print('= Examples =')
    print('Example 1 -- `dir` seqator mode.')
    print('  Move `.dna` files having contig coverage greater than 10 from directory `indir/` to `outdir/`.')
    print('  python3 seqator.py \\')
    print('    -i indir \\')
    print('    -x dna \\')
    print('    -p cov \\')
    print('    -f gt \\')
    print('    -t 10 \\')
    print('    -o outdir')
    print('Example 2 -- `fasta_file` seqator mode')
    print('  Copy sequences from file `input.fasta` which have sequence length less than 1000 bp to file `output.fasta.gz`.')
    print('  python3 seqator.py \\')
    print('    -i input.fasta \\')
    print('    -p len \\')
    print('    -f lt \\')
    print('    -t 1000 \\')
    print('    -o output.fasta.gz')
 
# end def

# Firstly check for information-providing flags

if '-h' in sys.argv[1:] or '--help' in sys.argv[1:]:
    print_help()
    platf_depend_exit()
# end if

if '-v' in sys.argv[1:] or '--version' in sys.argv[1:]:
    print(__version__)
    platf_depend_exit()
# end if


import argparse

parser = argparse.ArgumentParser()

# Input

parser.add_argument(
    '-i',
    '--input',
    help='TODO: add help',
    required=True
)

# Parameters

parser.add_argument(
    '-m',
    '--seqator-mode',
    help='TODO: add help',
    required=False,
    default='auto'
)

parser.add_argument(
    '-p',
    '--filter-parameter',
    help='TODO: add help',
    required=False,
    default='cov'
)

parser.add_argument(
    '-f',
    '--filter-mode',
    help='TODO: add help',
    required=False,
    default='lt'
)

parser.add_argument(
    '-t',
    '--threshold',
    help='TODO: add help',
    required=True
)

parser.add_argument(
    '-x',
    '--target-file-extention',
    help='TODO: add help',
    required=False,
    default='dna'
)

# Output

parser.add_argument(
    '-o',
    '--output',
    help='TODO: add help',
    required=False,
    default=''
)

parser_args = parser.parse_args()


# == Import them now ==

import os
import re
import gzip
import shutil
import operator


# == Data ==

# SPAdes-like string
SPADES_PATTERN = re.compile(
    r'NODE_[0-9]+_length_[0-9]+_cov_[0-9\.\,]+'
)


class Args:
    def __init__(self, parser_args):
        self.seqator_mode = parse_seqator_mode(parser_args)
        self.input = parse_input(parser_args, self.seqator_mode)
        self.threshold = parse_threshold(parser_args)
        self.target_file_extention = r'{}'.format(parser_args.target_file_extention)
        self.filter_parameter = parse_parameter(parser_args)
        self.filter_mode = parse_filter_mode(parser_args)
        self.output_path = parse_output_path(parser_args, self.seqator_mode)
    # end def
# end class


class FastaRecord:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
    # end def
# end class


class NoTargerFilesError(Exception):
    pass
# end class


# == Functions ==

def parse_threshold(parser_args):
    try:
        threshold = float(
            parser_args.threshold.replace(',', '.') # specially for evil russian KGB
        )
    except ValueError:
        print('Error: threshold value is invalid: `{}`'.format(parser_args.threshold))
        print('It must be a number.')
        platf_depend_exit(1)
    # end try
    return threshold
# end def

def parse_parameter(parser_args):
    allowed_parameters = {'len', 'cov'}
    filter_parameter = parser_args.filter_parameter.lower()
    if not filter_parameter in allowed_parameters:
        print('Error: target parameter is invalid: `{}`'.format(filter_parameter))
        print('It must be one of the following: {}.'.format(allowed_parameters))
        platf_depend_exit(1)
    # end if
    return filter_parameter
# end def

def parse_filter_mode(parser_args):
    allowed_modes = {
        'lt', 'le',
        'gt', 'ge',
        'eq', 'ne',
    }
    filter_mode = parser_args.filter_mode.lower()
    if not filter_mode in allowed_modes:
        print('Error: comparison mode is invalid: `{}`'.format(filter_mode))
        print('It must be one of the following: {}.'.format(allowed_modes))
        platf_depend_exit(1)
    # end if
    return filter_mode
# end def

def parse_seqator_mode(parser_args):
    allowed_modes = {
        'dir',
        'fasta_file',
        'auto',
    }
    seqator_mode = parser_args.seqator_mode.lower()
    if not seqator_mode in allowed_modes:
        print('Error: comparison mode is invalid: `{}`'.format(seqator_mode))
        print('It must be one of the following: {}.'.format(allowed_modes))
        platf_depend_exit(1)
    # end if

    input_abspath = os.path.abspath(parser_args.input)

    if seqator_mode == 'auto':
        if os.path.isdir(input_abspath):
            seqator_mode = 'dir'
        elif os.path.isfile(input_abspath):
            seqator_mode = 'fasta_file'
        else:
            print('Error: cannot recognize input type: `{}`'.format(parser_args.input))
            platf_depend_exit(1)
        # end if
    # end if

    return seqator_mode
# end def

def parse_input(parser_args, seqator_mode):
    input_abspath = os.path.abspath(parser_args.input)

    if seqator_mode == 'fasta_file':
        if not os.path.isfile(input_abspath):
            print('Error: `{}` is not a file'.format(input_abspath))
            platf_depend_exit(1)
        # end if
    elif seqator_mode == 'dir':
        if not os.path.isdir(input_abspath):
            print('Error: `{}` is not a directory'.format(input_abspath))
            platf_depend_exit(1)
        # end if
    else:
        print('Error: invalid seqator mode: `{}`'.format(seqator_mode))
        platf_depend_exit(1)
    # end if

    return input_abspath
# end def

def parse_output_path(parser_args, seqator_mode):
    if parser_args.output != '':
        return os.path.abspath(parser_args.output)
    else:
        if seqator_mode == 'fasta_file':
            return os.path.join(
                os.getcwd(),
                '{}_{}_{}.fasta'.format(
                    parser_args.filter_parameter,
                    parser_args.filter_mode,
                    parser_args.threshold
                )
            )
        elif seqator_mode == 'dir':
            return os.path.join(
                os.getcwd(),
                '{}_{}_{}'.format(
                    parser_args.filter_parameter,
                    parser_args.filter_mode,
                    parser_args.threshold
                )
            )
        else:
            print('Error: invalid seqator mode: `{}`'.format(seqator_mode))
            platf_depend_exit(1)
        # end if
    # end if
# end def


def report_args(args):
    print('-- Arguments --')
    if args.seqator_mode == 'dir':
        print('Input:')
        print('  Input directory: {}'.format(args.input))
        print('Output:')
        print('  Output directory: {}'.format(args.output_path))
    elif args.seqator_mode == 'fasta_file':
        print('Input:')
        print('  Input file: {}'.format(args.input))
        print('Output:')
        print('  Output file: {}'.format(args.output_path))
    # end if
    print('Parameters:')
    print('  Seqator mode: {}'.format(args.seqator_mode))
    print('  Target parameter: {}'.format(args.filter_parameter))
    print('  Comparison mode: {}'.format(args.filter_mode))
    print('  Threshold: {}'.format(args.threshold))
    if args.seqator_mode == 'dir':
        print('  Target file extention: {}'.format(args.target_file_extention))
    elif args.seqator_mode == 'fasta_file':
        print('  Target file extention: not applicable')
    # end if

    print('-' * 15 + '\n')
# end def


def get_input_fpaths(args):
    input_dirpath = args.input

    spades_like_fpaths = filter(
        fpath_is_spades_like,
        os.listdir(input_dirpath)
    )

    fpaths_with_target_extention = filter(
        lambda f: f.endswith('.{}'.format(args.target_file_extention)),
        spades_like_fpaths
    )

    make_abspath = lambda basename: os.path.abspath(
        os.path.join(input_dirpath, basename)
    )
    return tuple(
        map(
            make_abspath,
            fpaths_with_target_extention
        )
    )
# end def

def fpath_is_spades_like(fpath):
    global SPADES_PATTERN
    return not SPADES_PATTERN.search(fpath) is None
# end def


def fpath_to_spades_like_str(fpath, target_file_extention):
    return os.path.basename(
        fpath
    ).replace(
        '.{}'.format(target_file_extention),
        ''
    )
# end def


def choose_parse_target_param_function(filter_parameter):
    if filter_parameter == 'cov':
        return parse_coverage
    else:
        return parse_len
    # end def
# end def

def parse_coverage(spades_like_str):
    check_spades_str_format(spades_like_str)
    try:
        cov_str = spades_like_str.split('_')[5].replace(',', '.')
        cov = float(cov_str)
    except ValueError:
        print('Error in sequence `{}`.'.format(spades_like_str))
        print('  Coverage is invalid: `{}`.'.format(cov_str))
        print('It must be a number.')
        platf_depend_exit(1)
    # end try
    return cov
# end def

def parse_len(spades_like_str):
    check_spades_str_format(spades_like_str)
    try:
        len_str = spades_like_str.split('_')[3]
        length = int(len_str)
    except ValueError:
        print('Error in sequence `{}`.'.format(spades_like_str))
        print('  Length is invalid: `{}`.'.format(len_str))
        print('It must be an integer number.')
        platf_depend_exit(1)
    # end try
    return length
# end def

def check_spades_str_format(string):
    global SPADES_PATTERN
    ok = not SPADES_PATTERN.search(string) is None
    if not ok:
        print('Error: invalid format of a SPAdes-like string: `{}`'.format(string))
        print('Here is an example of a valid SPAdes-like string:')
        print('  NODE_1_length_40000_cov_22.55')
        platf_depend_exit(1)
    # end def
# end def


def choose_comparison_function(filter_mode):
    if filter_mode == 'lt':
        return operator.lt
    elif filter_mode == 'le':
        return operator.le
    elif filter_mode == 'gt':
        return operator.gt
    elif filter_mode == 'ge':
        return operator.ge
    elif filter_mode == 'eq':
        return operator.eq
    elif filter_mode == 'ne':
        return operator.ne
    else:
        raise ValueError('Invalid mode: `{}`'.format(filter_mode))
    # end if
# end def


def fasta_records(fpath):
    first_seq = True
    line = None
    name, seq = None, ''

    if fpath.endswith('.gz'):
        open_this_fasta = gzip.open
    else:
        open_this_fasta = open
    # end if

    with open_this_fasta(fpath, 'rt') as infile:
        while line != '':
            line = infile.readline()

            if line.startswith('>'):
                if not first_seq:
                    yield FastaRecord(name, seq)
                    name, seq = None, ''
                else:
                    first_seq = False
                # end if

                name = line.strip()[1:]

            elif line.strip() != '':
                seq += line.strip()

            else:
                yield FastaRecord(name, seq)
                return
            # end if
        # end while
    # end with
# end def


def write_fasta_record(outfile, record):
    outfile.write(
        '>{}\n{}\n'.format(
            record.name,
            record.seq
        )
    )
# end def


def make_outdir(args):
    if args.seqator_mode == 'dir':
        if not os.path.isdir(args.output_path):
            make_directory(args.output_path)
        # end if
    elif args.seqator_mode == 'fasta_file':
        outdir_path = os.path.dirname(args.output_path)
        if not os.path.isdir(outdir_path):
            make_directory(outdir_path)
        # end if
    else:
        print('Error: invalid seqator mode: `{}`'.format(args.seqator_mode))
        platf_depend_exit(1)
    # end if
# end def

def make_directory(dirpath):
    try:
        os.makedirs(dirpath)
    except OSError as err:
        print('Error: cannot create directory `{}`'.format(outdir_path))
        print(err)
        platf_depend_exit(1)
    # end try
# end def


def remove_outdir(outdir_path):
    try:
        os.rmdir(outdir_path)
    except OSError as err:
        print('Error: cannot remove output directory `{}`'.format(outdir_path))
        print(err)
        platf_depend_exit(1)
    # end try
# end def


def bin_files_dir_mode(args):
    input_fpaths = get_input_fpaths(args)

    if len(input_fpaths) == 0:
        print("Warning: not appropriate input files found.")
        print('Advice: please, check target extention (option -x/--target-file-extention)')
        return 0, 0
    # end if

    parse_filter_parameter = choose_parse_target_param_function(args.filter_parameter)
    check_filter_parameter = choose_comparison_function(args.filter_mode)

    moved_counter = 0
    for fpath in input_fpaths:
        # Retrieve target parameter from file name
        spades_like_str = fpath_to_spades_like_str(fpath, args.target_file_extention)
        targ_param_value = parse_filter_parameter(spades_like_str)
        # Compare and move if matches
        if check_filter_parameter(targ_param_value, args.threshold):
            print('Moving `{}`'.format(os.path.basename(fpath)))
            shutil.move(fpath, args.output_path)
            moved_counter += 1
        # end if
    # end for

    # Remove output dir if it is empty
    if len(os.listdir(args.output_path)) == 0:
        print('Output directory is empty -- removing it')
        remove_outdir(args.output_path)
    # end if

    total_count = len(input_fpaths)
    return total_count, moved_counter
# end def


def bin_file_fasta_file_mode(args):
    parse_filter_parameter = choose_parse_target_param_function(args.filter_parameter)
    check_filter_parameter = choose_comparison_function(args.filter_mode)

    fasta_outfpath = args.output_path
    if fasta_outfpath.endswith('.gz'):
        open_this_fasta = gzip.open
    else:
        open_this_fasta = open
    # end if

    total_counter, cp_counter = 0, 0
    with open_this_fasta(fasta_outfpath, 'wt') as outfile:

        for record in fasta_records(args.input):
            total_counter += 1
            # Retrieve target parameter from file name
            targ_param_value = parse_filter_parameter(record.name)
            # Compare and copy if matches
            if check_filter_parameter(targ_param_value, args.threshold):
                print('Copying sequence `{}`'.format(record.name))
                write_fasta_record(outfile, record)
                cp_counter += 1
            # end if
        # end for
    # end with

    return total_counter, cp_counter
# end def


# == Proceed ==

print('\n  == seqator version {} == \n'.format(__version__))

# Parse and validate arguments
args = Args(parser_args)
report_args(args)

# Create output dir
make_outdir(args)

# Binning
print('Start binning')
if args.seqator_mode == 'dir':
    total_count, moved_count = bin_files_dir_mode(args)
    print('Binning is completed')
    print('=' * 45)
    print('{}/{} files moved to `{}/`'.format(moved_count, total_count, args.output_path))
    print('{} files left in `{}/`'.format(total_count - moved_count, args.input))
elif args.seqator_mode == 'fasta_file':
    total_count, cp_count = bin_file_fasta_file_mode(args)
    print('Binning is completed')
    print('=' * 45)
    print(
        '{}/{} sequences copied from file `{}` to file `{}`' \
            .format(cp_count, total_count, args.input, args.output_path)
    )
# end if

print('\nCompleted!')
platf_depend_exit(0)
