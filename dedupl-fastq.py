#!/usr/bin/env python3
# -*- cofing: utf-8 -*-

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
    print("\nPython {}.{}.{}".format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro))
# end if

__version__ = "1.0.a"
# Last updated: 29.06.2020


import os

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


outdir = os.path.join(os.getcwd(), "fastq_deduplicated")

print("  == dedupl_fastq.py (Version {}) ==".format(__version__))
print("\nThis script deduplicates fastq files in the working directory.")
print("Deduplicated files will be written to directory '{}'".format(outdir))

if len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help", "-help"):
    print("\n Usage:")
    print("No arguments, just run it to process all fastq files in the working directory:")
    print("  ./dedupl-fastq.py")
    platf_depend_exit(0)
# end if

import glob
from re import search as re_search
from gzip import open as open_as_gzip


from time import time, strftime, gmtime

start_time = time() # consider time of importing as start time

def getwt():
    return strftime("%H:%M:%S", gmtime( time() - start_time))
# end def getwt


is_fastq = lambda f: True if not re_search(r".+\.f(ast)?q(\.gz)?$", f) is None else False
infiles = sorted(tuple(filter(is_fastq, os.listdir('.'))))

if len(infiles) == 0:
    print("\nNo fastq files in the working directory")
    print("\n Usage:")
    print("No arguments, just run it to process all fastq files in the working directory:")
    print("  ./dedupl-fastq.py")
    if sys.platform.startswith("win"):
        input("Press ENTER to exit")
    # end if
    exit(1)
else:
    print("\n{} files are found and will be processed.".format(len(infiles)))
# end if


id_lst = list()
glob_uniq, glob_dupl = 0, 0

write_mode = 'w'

if not os.path.isdir(outdir):
    try:
        os.makedirs(outdir)
    except OSError as oserr:
        print("Cannot create output directory: {}".format( str(oserr) ))
        platf_depend_exit(1)
    # end try
# end if

if len(os.listdir(outdir)) != 0:
    print("\nOutput directory is not empty!")
    error = True
    while error:
        reply = input("""  Overwrite/append/quite? [O/a/q] """)
        if reply.upper() == "O":
            error = False
            for fpath in glob.iglob(os.path.join(outdir, '*')):
                print("Removing '{}'.".format(fpath))
                try:
                    os.unlink(fpath)
                except OSError as oserr:
                    print("Cannot remove file '{}': {}".format( fpath, str(oserr) ))
                    platf_depend_exit(1)
                # end try
            # end for
            print()
        elif reply == 'a':
            error = False
            write_mode = 'a'
        elif reply == 'q':
            sys.exit(0)
        else:
            print("Invalid reply: '{}'".format(reply))
        # end if
    # end while
# end if


OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode("utf-8").strip()  # format gzipped line
)

is_gzipped = lambda file: True if file.endswith(".gz") else False

# GNU gzip utility is faster, but there can be presence of absence of it :)
gzip_util = "gzip"
util_found = False
for directory in os.environ["PATH"].split(os.pathsep):
    if os.path.isdir(directory) and gzip_util in os.listdir(directory):
        util_found = True
        break
    # end if
# end for

if util_found:

    def gzip_func(fpath):
        """Function that compresses output FASTA and FASTQ files with gzip utility.
        
        :param fpath: path to file to compress;
        :type fpath: str;
        """

        try:
            os.system("{} {}".format(gzip_util, fpath))
        except OSError as oserr:
            print("cannot gzip file '{}'".format(os.path.basename(fpath)))
            print("Reason: {}".format( str(oserr) ))
            # Try to gzip others -- continue loop
        else:
            print("{} - Gzipping is completed".format(getwt()))
        # end try
    # end def gzip_func
else:

    from shutil import copyfileobj as shutil_copyfileobj

    def gzip_func(fpath):
        """Function that compresses output FASTA and FASTQ files with Python gzip and shutil modules.
        
        :param fpath: path to file to compress;
        :type fpath: str;
        """
        try:
            # form .fasta.gz file 'by hand'
            with open(fpath, 'rb') as fastqa_file, open_as_gzip(fpath+".gz", "wb") as faqgz_file:
                shutil_copyfileobj(fastqa_file, faqgz_file)
            # end with
            os.unlink(fpath) # remove plain FASTA file
        except OSError as oserr:
            print("cannot gzip file '{}'".format(os.path.basename(fpath)))
            print("Reason: {}".format( str(oserr) ))
            # Try to gzip others -- continue loop
        else:
            print("{} - Gzipping is completed".format(getwt()))
        # end try
    # end def gzip_func
# end if


print("{} - Start deduplication\n".format(getwt()))


for fq_path in infiles:

    print("{} - Processing '{}'...".format(getwt(), fq_path))

    uniq, dupl = 0, 0

    new_path = os.path.basename(fq_path.replace(".fastq", "_dedupl.fastq"))
    new_path = new_path.replace(".gz", "")
    new_path = os.path.join(outdir, new_path)

    how_to_open = OPEN_FUNCS[ is_gzipped(fq_path) ]
    fmt_func = FORMATTING_FUNCS[ is_gzipped(fq_path) ]

    with how_to_open(fq_path) as infile, open(new_path, write_mode) as outfile:

        while True:

            seq_id = fmt_func(infile.readline())

            if seq_id == "":
                break
            # end if

            if seq_id not in id_lst:
                id_lst.append(seq_id)

                outfile.write(seq_id + '\n')
                for _ in range(3):
                    outfile.write(fmt_func(infile.readline()) + '\n')
                # end for
                uniq += 1
            else:
                for _ in range(3):
                    infile.readline()
                # end for
                dupl += 1
            # end if
        # end while
    # end with

    glob_uniq += uniq
    glob_dupl += dupl

    print("{} unique sequences found.".format(uniq))
    if dupl != 0:
        print("{} duplicated sequences found and removed.".format(dupl))
    # end if

    if uniq != 0:
        if is_gzipped(fq_path):
            print("{} - Gzipping result file '{}'...".format(getwt(), os.path.basename(new_path)))
            gzip_func(new_path)
            new_path += ".gz"
        # end if
    # end if
    print('-'*20)
# end for

print("\n{} - Deduplication is completed".format(getwt()))
print("Finally:")
print("  {} unique sequences found.".format(glob_uniq))
print("  {} duplicated sequences found and removed.\n".format(glob_dupl))

platf_depend_exit(0)