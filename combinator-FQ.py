#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

__version__ = "1.0.a"
# Year, month, day
__last_update_date__ = "2020-03-02"

# |===== Check python interpreter version =====|

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


def err_fmt(text=""):
    """Function for configuring error messages"""
    return "\n  \a!! - ERROR: " + text + '\n'
# end def err_fmt


if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    print("combinator-FQ. Version {}. {} edition.".format(__version__, __last_update_date__))
    print("Script screens ends of contigs for identity in order to facilitate further manual scaffolding.")
    print("\nIt generates output table 'combinator_output_FQ.tsv' of following format:")
    print("#  Contig name  Length  Coverage  GC(%)  Multiplicity  Annotation  Start      End")
    print("1  NODE_2...    304356  93.7155   34.12  1.2                       S_NODE_22  S_NODE_26")
    platf_depend_exit(0)
# end if

if "-v" in sys.argv[1:] or "--version" in sys.argv[1:]:
    print(__version__)
    platf_depend_exit(0)
# end if


import getopt

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "hvp:k:i:a:d:",
        ["help", "version", "prefix=", "k-mer-len=", "mink=", "maxk=", "indir="])
except getopt.GetoptError as gerr:
    print( str(gerr) )
    platf_depend_exit(2)
# end try

import re
import os
import glob

prefix = ""
mink = 21
maxk = 127

is_fasta = lambda f: not re.search(r"(m)?f(ast)?a(\.gz)?$", f) is None


if len(args) == 0:

    cont_fpaths = glob.glob( os.path.join(os.getcwd(), "*contigs.f*") )
    cont_fpaths = tuple(filter( is_fasta, cont_fpaths ))

    if len(cont_fpaths) != 0:
        contigs_fpath = cont_fpaths[0]
        print("File '{}' is found and will be processed.".format(contigs_fpath))
        error = True
        while error:
            reply = input("""Press ENTER to process this file
  or enter 'q' to quit:>> """)
            if reply == "":
                error = False
            elif reply == 'q':
                sys.exit(0)
            else:
                print("Invalid reply: '{}'\n".format(reply))
            # end if
        # end while

    else:
        print("No files considered as input fot combinator-FQ is are found.")
        print("Exitting...")
        platf_depend_exit(1)
    # end if

else:
    contigs_fpath = args[0]
    if not os.path.exists(contigs_fpath):
        print("File '{}' does not exist.".format(contigs_fpath))
        platf_depend_exit(1)
    # end if
# end if

for opt, arg in opts:

    if opt in ("-k", "--k-mer-len"):

        try:
            arg = int(arg)
            if arg <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\nLength of k-mer must be positive integer number.")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        else:
            mink, maxk = k, k
        # end try

    elif opt in ("-i", "--mink"):

        try:
            mink = int(arg)
            if mink <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\nMinimum length of k-mer must be positive integer number.")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-a", "--maxk"):

        try:
            maxk = int(arg)
            if maxk <= 0:
                raise ValueError
            # end if
        except ValueError:
            print("\nMaximum length of k-mer must be positive integer number.")
            print("Your value: '{}'".format(arg))
            platf_depend_exit(1)
        # end try

    elif opt in ("-p", "--prefix"):
        prefix = arg
    # end if
# end for


if mink > maxk:
    print(err_fmt("Minimum length of k-mer is greater than maximum length of k-mer."))
    print("Values specified by you:")
    print("Minimum length of k-mer: {}.".format(mink))
    print("Maximum length of k-mer: {}.".format(maxk))
    platf_depend_exit(1)
# end if


# Dictionary maps complementary bases according to IUPAC:
compl_dict = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
    'H': 'D', 'V': 'B', 'U': 'A', 'N': 'N',
    'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
    'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
    'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
    'h': 'd', 'v': 'b', 'u': 'a', 'n': 'n'
}

# Function returns complement "comrade" of given base
get_compl_base = lambda base: compl_dict[base]

# Function returns reverse-complement "comrade" of passed DNA sequence
rc = lambda seq: "".join(map( get_compl_base, reversed(seq) ))


from gzip import open as open_as_gzip

is_gzipped = lambda file: True if file.endswith(".gz") else False

# For opening plain text and gzipped files
OPEN_FUNCS = (open, open_as_gzip)

# Data from plain text and gzipped should be parsed in different way,
#   because data from .gz is read as 'bytes', not 'str'.
FORMATTING_FUNCS = (
    lambda line: line.strip(),   # format text line
    lambda line: line.decode("utf-8").strip()  # format gzipped line
)


print("\nFile '{}' if processing...\n".format(contigs_fpath))

open_func = OPEN_FUNCS[ is_gzipped(contigs_fpath) ]
fmt_func = FORMATTING_FUNCS[ is_gzipped(contigs_fpath) ]

contigs_names = list() # list for names of contigs
contigs_seqs = list() # list for sequences of contigs
seq = ""

with open_func(contigs_fpath) as contigs_file:

    for line in contigs_file: # line by line
        line = fmt_func(line)
        if line[0].startswith('>'):
            contigs_names.append(line)
            contigs_seqs.append(seq)
            seq = ""
        else:
            seq += line
        # end if
    # end for
# end with

contigs_seqs.append(seq) # append last sequence to list
contigs_seqs = contigs_seqs[1:] # remove first empty string

print(len(contigs_names))
print(len(contigs_seqs))

exit(0)

voc_ends = dict() # dictionary for collecting contigs data
# Variables for collecting statistics
total_length = 0
agv_coverage = 0
min_coverage = float("inf")
max_coverage = 0

# Readable indices
FULL_NAME, NAME, LEN, COV, GC, START, RC_START, END, RC_END, START_MATCH, END_MATCH = range(11)

# Open output file
combinator_output = open("{}_combinator_output_FQ.tsv".format(prefix), 'w')
# Write head of the table:
combinator_output.write("#\tContig name\tLength\tCoverage\tGC(%)\tMultiplicity\tAnnotation\tStart\tEnd\n")

# Iterate over contigs and form voc_ends dictionary
for i, contig_name in enumerate(contigs_names):

    print(contig_name)

    contig_len = len(contigs_seqs[i])

    # Calculating GC-content
    gc_content = 0
    for up_base, low_base in zip(('G', 'C', 'S'),('g', 'c', 's')):
        gc_content += contigs_seqs[i].count(up_base) + contigs_seqs[i].count(low_base)
    # end for

    gcContent = round((gc_content / len(contigs_seqs[i]) * 100), 2)

    # Parse fasta header:
    cov = str(contig_name.split('_')[5]) # get coverage
    name = 'NODE_' + contig_name.split('_')[1] # get name in 'NODE_<NUMBER>' format

    # Some common info
    voc_ends.setdefault(i, []).append(contig_name[1:])              # FULL_NAME
    voc_ends.setdefault(i, []).append(name)                         # NAME
    voc_ends.setdefault(i, []).append(contig_len)                   # LEN
    voc_ends.setdefault(i, []).append(cov)                          # COV
    voc_ends.setdefault(i, []).append(gcContent)                    # GC

    # Ends of contigs:
    voc_ends.setdefault(i, []).append(contigs_seqs[i][:maxk])       # START
    voc_ends.setdefault(i, []).append(rc(contigs_seqs[i][:maxk]))   # RC_START
    voc_ends.setdefault(i, []).append(contigs_seqs[i][-maxk:])      # END
    voc_ends.setdefault(i, []).append(rc(contigs_seqs[i][-maxk:]))  # EC_END

    # Ends' matches
    voc_ends.setdefault(i, []).append('-')                          # START_MATCH
    voc_ends.setdefault(i, []).append('-')                          # END_MATCH

    # Collecting statistics
    flt_cov = float(cov)

    total_length += contig_len
    agv_coverage += flt_cov

    min_coverage = min(min_coverage, flt_cov)
    max_coverage = max(max_coverage, flt_cov)
# end for

print('\n' + '=' * 45)

# поиск совпадений между концами контигов
print(u'Поиск совпадений между концами контигов...')
log = u'\n Совпадения между концами контигов:'
i3 = 0
while i3 < i:
    i32 = 0
    print('\n' + '-' * 45 + '\n' + voc_ends[i3][1] + ': ')
    log += ('\n' + '-' * 45 + '\n' + voc_ends[i3][1] + ': ')
    # Проверка комплементарности НАЧАЛА нода с другими концами
    while i32 < i:
        # 1-1
        if i32 != i3 and voc_ends[i3][5] == voc_ends[i32][5]:
            print(u' 1-1: Начало ' + voc_ends[i3][1] + u' совпадает c началом ' + voc_ends[i32][1])
            log += (u'\n 1-1: Начало ' + voc_ends[i3][1] + u' совпадает c началом ' + voc_ends[i32][1])
        # 1-2
        if voc_ends[i3][5] == voc_ends[i32][6]:
            if voc_ends[i3][9] == '-':
                voc_ends[i3][9] = ('H_' + voc_ends[i32][1])
            else:
                voc_ends[i3][9] += ('  H_' + voc_ends[i32][1])
            print(u' 1-2: Начало ' + voc_ends[i3][1] + u' совпадает c ОК-началом ' + voc_ends[i32][1])
            log += (u'\n 1-2: Начало ' + voc_ends[i3][1] + u' совпадает c ОК-началом ' + voc_ends[i32][1])
        # 1-3
        if voc_ends[i3][5] == voc_ends[i32][7]:
            if i3 == i32:
                if voc_ends[i3][9] == '-':
                    voc_ends[i3][9] = u'|кольцо|'
                else:
                    voc_ends[i3][9] += u'|замыкается в кольцо|'
                print(' 1-3: ' + voc_ends[i3][1] + u' замыкается в кольцо!')
                log += ('\n 1-3: ' + voc_ends[i3][1] + u' замыкается в кольцо!')
            else:
                if voc_ends[i3][9] == '-':
                    voc_ends[i3][9] = (voc_ends[i32][1] + '_K')
                else:
                    voc_ends[i3][9] += ('  ' + voc_ends[i32][1] + '_K')
                print(u' 1-3: Начало ' + voc_ends[i3][1] + u' совпадает c концом ' + voc_ends[i32][1])
                log += (u'\n 1-3: Начало ' + voc_ends[i3][1] + u' совпадает c концом ' + voc_ends[i32][1])
        # 1-4
        if voc_ends[i3][5] == voc_ends[i32][8]:
            if i3 == i32:
                print(u' 1-4: Начало ' + voc_ends[i3][1] + u' идентично его RC концу!')
                log += (u'\n 1-4: Начало ' + voc_ends[i3][1] + u' идентично его RC концу!')
            else:
                print(u' 1-4: Начало ' + voc_ends[i3][1] + u' совпадает c ОК-концом ' + voc_ends[i32][1])
                log += (u'\n 1-4: Начало ' + voc_ends[i3][1] + u' совпадает c ОК-концом ' + voc_ends[i32][1])
        i32 += 1

    # Проверка комплементарности КОНЦА нода с другими концами
    i32 = 0
    while i32 < i:
        # 3-1
        if voc_ends[i3][7] == voc_ends[i32][5] and i3 != i32:
            if voc_ends[i3][10] == '-':
                voc_ends[i3][10] = ('H_' + voc_ends[i32][1])
            else:
                voc_ends[i3][10] += ('  H_' + voc_ends[i32][1])
            print(u' 3-1: Конец ' + voc_ends[i3][1] + u' совпадает c началом ' + voc_ends[i32][1])
            log += (u'\n 3-1: Конец ' + voc_ends[i3][1] + u' совпадает c началом ' + voc_ends[i32][1])
        # 3-2
        if i3 != i32 and voc_ends[i3][7] == voc_ends[i32][6]:
            print(u' 3-2: Конец ' + voc_ends[i3][1] + u' совпадает c ОК-началом ' + voc_ends[i32][1])
            log += (u'\n 3-2: Конец ' + voc_ends[i3][1] + u' совпадает c ОК-началом ' + voc_ends[i32][1])
        # 3-3
        if i32 != i3 and voc_ends[i3][7] == voc_ends[i32][7]:
            print(u' 3-3: Конец ' + voc_ends[i3][1] + u' совпадает c концом ' + voc_ends[i32][1])
            log += (u'\n 3-3: Конец ' + voc_ends[i3][1] + u' совпадает c концом ' + voc_ends[i32][1])
        # 3-4
        if voc_ends[i3][7] == voc_ends[i32][8]:
            if voc_ends[i3][10] == '-':
                voc_ends[i3][10] = (voc_ends[i32][1] + '_K')
            else:
                voc_ends[i3][10] += ('  ' + voc_ends[i32][1] + '_K')
            print(u' 3-4: Конец ' + voc_ends[i3][1] + u' совпадает c ОК-концом ' + voc_ends[i32][1])
            log += (u'\n 3-4: Конец ' + voc_ends[i3][1] + u' совпадает c ОК-концом ' + voc_ends[i32][1])
        i32 += 1

    i3 += 1
# Формирование итоговой таблицы в файле combinator_output
i3 = 0
qqq = 0 # переменная для рассчёта LQ-коэфициэнта
while i3 < i:
    combinator_output.write(str(i3 + 1) + '\t')
    combinator_output.write(voc_ends[i3][0] + '\t')
    combinator_output.write(str(voc_ends[i3][2]) + '\t')
    combinator_output.write(str(voc_ends[i3][3].replace('.', ',')) + '\t')
    combinator_output.write(voc_ends[i3][4] + '\t')
    combinator_output.write(str(round(float(voc_ends[i3][3])/float(voc_ends[0][3]), 1)).replace('.', ',') + '\t\t')  # расчёт коэффицента кратности контига
    combinator_output.write(voc_ends[i3][9] + '\t' + voc_ends[i3][10] + '\n')
    # сбор данных для рассчёта LQ-коэфициэнта
    if voc_ends[i3][9] == '-':
        qqq += 1
    elif voc_ends[i3][9] == '|кольцо|':
        qqq -= 1
    if voc_ends[i3][10] == '-':
        qqq += 1
    i3 += 1

# запись статистической информации в файл
combinator_output.write('\n' + '=' * 45 + '\n')
combinator_output.write('Total Length:   \t' + str(total_length) + '\n')
combinator_output.write('Min Coverage:   \t' + str(min_coverage) + '\n')
combinator_output.write('Max Coverage:   \t' + str(max_coverage) + '\n')
combinator_output.write('Mid Coverage:   \t' + str(round(agv_coverage / i, 3)) + '\n')
combinator_output.write('Total contigs processed: \t' + str(i) + '\n')
#combinator_output.write('QQQ:   \t' + str(qqq) + '\n')
combinator_output.write('LQ-coefficient: \t' + str((round(100-(qqq*100)/(i*2), 2))) + '\n')
combinator_output.write('=' * 45)

combinator_output.write('\n' + log)     # дописывание лога в конец файла.
combinator_output.close()               # закрываем файл
print('\n' + '=' * 45)
print(u'Готово!')
input(u'Для завершения программы, нажмите \"ENTER\"...')
exit()