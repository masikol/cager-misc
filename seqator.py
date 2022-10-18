#!/usr/bin/env python3
#
# ---------------------------------
# Script moves "SPAdes-like" *.dna files with coverage less than specified one
#   from 'contigs/' directory to directory 'cov_below_x/'.
# "SPAdes-like" means that name of file is of following format:
#   NODE_1_length_61704_cov_114.517.dna
# If directory `./contigs/` contains no "SPAdes-like" files,
#   seqator.py takes them from directory `./contigs/DNA Files`
# ---------------------------------

__version__ = "2.0.a"
__last_update_date__ = "2022-10-18"


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
# end if


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
# end def

def print_help():
    print("\n    |=== seqator ===|")
    print("Version {}. {} edition.\n".format(__version__, __last_update_date__))

    print("Script moves \"SPAdes-like\" *.dna files with coverage less than specified one")
    print("  from directory `./contigs/` to directory `cov_below_x/`.")
    print("\"SPAdes-like\" means that name of file is of following format:")
    print("  NODE_1_length_61704_cov_114.517.dna")
    print("If directory `./contigs/` contains no \"SPAdes-like\" files,")
    print("  seqator.py takes them from directory `./contigs/DNA Files`.")

    print("\nMoreover, seqator performs the same task on file `./contigs.fasta`,")
    print("  if it is located in the working directory.")
    print("The script copies sequences having coverage below the threshold to a separate fasta file.")

    print("\nUsage:")
    print("You can specify coverage threshold in command line:\n")
    print("  python3 seqator.py 12\n")
    print("Or run just this, and script will ask you to enter threshold from keyboard:")
    print("  python3 seqator.py")
# end def

# Firstly check for information-providing flags

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print_help()
    platf_depend_exit()
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit()
# end if

# Get coverage value
if len(sys.argv) == 2:
    # If it is specified in argv

    try:
        cov_threshold = float(sys.argv[1].replace(',', '.')) # specially for evil russian KGB
        if cov_threshold < 0:
            raise ValueError
    except ValueError:
        print("Invalid coverage specified: `{}`".format(cov_threshold))
        print("Must be positive number!")
        platf_depend_exit(1)
    # end try

elif len(sys.argv) == 1:
    # If no coverage specified in argv, ask to enter it from keyboard

    cov_threshold = None
    while cov_threshold is None:
        cov_threshold = input("Please, enter coverage threshold:")
        cov_threshold = cov_threshold.replace(',', '.') # specially for evil russian KGB
        try:
            cov_threshold = float(cov_threshold)
            if cov_threshold < 0:
                raise ValueError
        except ValueError:
            print("Invalid coverage specified: `{}`".format(cov_threshold))
            print("It must be a positive number.")
            cov_threshold = None
        # end try
    # end while

else:
    print_help()
    platf_depend_exit(1)
# end if


import os
import re
import shutil


# We will process only those .dna files, which are SPAdes contigs
is_target_file = lambda f: not re.match(r"NODE_[0-9]+_length_[0-9]+_cov_[0-9\.\,]+\.dna", f) is None


def get_target_dirpath():
    possible_target_dirs = (
        "contigs",
        os.path.join("contigs", "DNA Files")
    )
    target_file_counter = 0
    target_dirpath = possible_target_dirs[target_file_counter]

    target_file_count = count_target_files(target_dirpath)

    while target_file_count == 0:
        target_file_counter += 1
        try:
            target_dirpath = possible_target_dirs[target_file_counter]
        except IndexError:
            raise NoDNAFilesError()
        # end try
        target_file_count = count_target_files(target_dirpath)
    # end while

    return target_dirpath
# end def


def count_target_files(target_dirpath):
    return len(
        tuple(
            filter(is_target_file, os.listdir(target_dirpath))
        )
    )
# end def


class NoDNAFilesError(Exception):
    pass
# end class


def get_target_files(target_dirpath):

    make_abspath = lambda basename: os.path.abspath(
        os.path.join(target_dirpath, basename)
    )

    return tuple(
        map(
            make_abspath,
            filter(
                is_target_file,
                os.listdir(target_dirpath)
            )
        )
    )
# end def


def parse_coverage(spades_like_str):
    # Get coverage + (maybe) `.dna`
    cov = spades_like_str.split("_")[5].replace(",", ".")
    # Prune '.dna'
    cov = cov.replace(".dna", "")
    return cov
# end def


def fasta_records(fpath):

    first_seq = True
    line = None
    name, seq = None, ""

    with open(fpath, "rt") as infile:

        while line != "":
            line = infile.readline()

            if line.startswith(">"):
                if not first_seq:
                    yield {
                        "name": name,
                        "seq" :  seq,
                    }
                    name, seq = None, ''
                else:
                    first_seq = False
                # end if

                name = line.strip()[1:]

            elif line.strip() != '':
                seq += line.strip()

            else:
                yield {
                    "name": name,
                    "seq" :  seq,
                }
                return
            # end if
        # end while
    # end with
# end def


def write_fasta_record(outfile, record):
    outfile.write('>{}\n{}\n'.format(
        record["name"],
        record["seq"])
    )
# end def



# Get paths to files to process
try:
    target_dirpath = get_target_dirpath()
except NoDNAFilesError:
    input_fpaths = tuple()
else:
    # take only 'NODE_*.dna' files
    input_fpaths = get_target_files(target_dirpath)
# end try

contigs_fpath = "contigs.fasta"
contigs_exist = os.path.exists(contigs_fpath)


if len(input_fpaths) == 0 and not contigs_exist:
    print("No input data found!")
    print('Cannot find appropriate ("SPAdes-like") .dna files,')
    print("  and the file `{}` does not exist.".format(contigs_fpath))
    print("Please, type `python3 seqator.py -h` to see help.")
    platf_depend_exit(1)
# end if


print(" seqator.py; Version {}; {} edition;\n".format(__version__, __last_update_date__))


# == Move DNA files ==

if len(input_fpaths) != 0:
    print("seqator will process directory `{}`".format(target_dirpath))

    # Configure path to output dir
    new_path = os.path.join(target_dirpath, "cov_below_{}".format(cov_threshold))
    if not os.path.isdir(new_path):
        os.makedirs(new_path)
    # end if
    print("Sequences with coverage below {} will be moved to directory `{}`\n".format(cov_threshold, new_path))

    # Proceed
    moved_counter = 0
    for fpath in input_fpaths:

        # Retrieve coverage from file name
        cov = parse_coverage(os.path.basename(fpath))

        # Compare and move if less
        if float(cov) < cov_threshold:
            print("Moving `{}`".format(fpath))
            shutil.move(fpath, new_path)
            moved_counter += 1
        # end if
    # end for


    # Remove output dir if it is empty
    if len(os.listdir(new_path)) == 0:
        os.rmdir(new_path)
    # end if

    print("=" * 45)
    print("{} files moved to `{}/`".format(moved_counter, new_path))
    print("{} files left in `{}/`\n".format(len(input_fpaths) - moved_counter, target_dirpath))
else:
    print("No .dna files to process have been found.\n")
# end if


# == Perform binning of file `contigs.fasta` if it exists in the WD ==

if contigs_exist:
    fasta_outfpath = os.path.join(os.getcwd(), "contigs_cov_below_{}.fasta".format(cov_threshold))
    print(
        "Copying sequences with coverage < {} from file `{}` to file `{}`." \
            .format(cov_threshold, contigs_fpath, fasta_outfpath)
    )

    total_counter, cp_counter = 0, 0
    with open(fasta_outfpath, "wt") as outfile:
        for record in fasta_records(contigs_fpath):
            cov = parse_coverage(record["name"])
            total_counter += 1

            if float(cov) < cov_threshold:
                print("Copying sequence `{}`".format(record["name"]))
                write_fasta_record(outfile, record)
                cp_counter += 1
            # end if
        # end for
    # end with

    print("=" * 45)
    print(
        "{}/{} sequences copied from file `{}` to file `{}`\n" \
            .format(cp_counter, total_counter, contigs_fpath, fasta_outfpath)
    )
else:
    print("No `contigs.fasta` file has been found.\n")
# end if


print("Completed!")
platf_depend_exit(0)
