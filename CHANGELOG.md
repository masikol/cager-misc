# cager-misc changelog

## 2025-03-12

Remove useless `genome_id_file_name_map` from `AAI-taxon.sh`.

### Version changes:

- AAI-taxon.sh: `1.0.a -> 1.0.b`

## 2025-03-07. Evening edition

Add `AAI-taxon.sh` v1.0.a.

## 2025-03-07

Add `gla-glar.py` v1.0.a.

## 2024-12-31

Add scripts `cigar_maplen.py` v1.0.a and `samtool_setop.py` v1.0.a.

## 2024-03-27

Add verbose mode to `sum-up-snv.awk`

### Version changes:

- sum-up-snv.awk: `0.1.c -> 0.2.a`

## 2023-03-01

Rewrite `seqator.py` compeletely to add flexibility.

### Version changes:

- seqator: `2.0.a -> 3.0.a`

## 2021-08-10

`most-freq-subseq` can now count subsequences either on positive strand or on both strands of input sequences (see option `--both-strands`). Default behaviour is to count subsequences only on positive strand.

### Version changes:

- most-freq-subseq: `1.1.b -> 1.2.a`

## 2021-04-16

seqator now takes input files from directory `./contigs/DNA Files` if directory `./contigs` contains no "SPAdes-like" `.dna` files.

### Version changes:

- seqator: `1.0.a -> 1.1.a`

## 2021-03-26

Moved combinator-FQ to the separate repository: [https://github.com/masikol/combinator-FQ](https://github.com/masikol/combinator-FQ).

## 2021-01-21 edition

- fasta-GC-content: output files now contains "Coverage" column.
- fasta-GC-content: renames column "S (G or C)" to "S (degenerate)".

### Version changes:

- fasta-GC-content: `1.1.a -> 1.1.b -> 1.1.c`

## 2021-01-21 edition

- combinator-FQ: bug fix (combinator used to ternimate in the end if it couldn't extract ordinal number of contig from contig's name);

### Version changes:

- combinator-FQ: `1.3.e -> 1.3.f`

## 2021-01-06 edition

- kromsatel: fixed bug that would cause kromsatel to prefer minor amplicons over major ones.

### Version changes:

- kromsatel: `1.2.e -> 1.3.a`

## 2020-12-03 edition

- combinator-FQ: `-o` option addded.

### Version changes:

- combinator-FQ: `1.3.d -> 1.3.e`

## 2020-12-02 edition

- kromsatel: output files naming changed: now "cleaned" is suffix, not prefix.

### Version changes:

- kromsatel: `1.2.d -> 1.2.e`

## 2020-11-14 edition

- kromsatel now works 2 times faster.

### Version changes:

- kromsatel: `1.2.c -> 1.2.d`

## 2020-11-13 edition

- kromsatel: performance improved.

### Version changes:

- kromsatel: `1.2.b -> 1.2.c`

## 2020-11-12 edition

- kromsatel version `1.2.a` added, then added version `1.2.b` with code commented and some optimized operations.

## 2020-10-23 edition

- pub: checked algorithm of finding color ballanced set of barcodes: now pub reclusters all barcodes except of the one having minimal silhouette score.

### Version changes:

- pub: `1.0.a --> 1.1.a`

## 2020-10-15 edition

- sum-up-snv (version `0.1.c`) added.

## 2020-10-08 edition

- combinator-FQ now rounds coverage with 2 trailing digits.

### Version changes:

- combinator-FQ: `1.3.c --> 1.3.d`

## 2020-06-29 edition

- dedupl-fastq.py added.

## 2020-06-14 edition

- pub.R added.

## 2020-06-15 edition.

- combinator-FQ: errorneous handling of SPADEs's `NODE_1` with zero coverage fixed;

### Version changes:

- combinator-FQ: `1.3.b --> 1.3.c`

## 2020-04-24 edition.

- combinator-FQ: compatibility with a5 improved;

### Version changes:

- combinator-FQ: `1.3.a --> 1.3.b`

## 2020-04-16 edition.

- fastq-read-count: now counts number of bases additionaly;

### Version changes:

- fastq-read-count: `1.0.a --> 1.1.a`;

## 2020-04-15 edition.

- combinator-FQ: calculating of expected genome length is improved;

### Version changes:

- combinator-FQ: `1.2.b --> 1.3.a`;

## 2020-04-14 edition.

- find-seq script added;
- NOS added;
- packer-dna-to-fasta added;
- seqator added;

## 2020-04-13 edition.

- fasta-GC-content: format of output file changed to tab-separated table with summary in the end;
- mean-qual script added;
- fastq2fasta script added;
- fastq-read-count script added;
- dna-summary script added;
- most-freq-subseq script added;

## 2020-04-10 evening edition.

- combinator-FQ: a5-compatibility bug fixed;

### Version changes:

1. combinator-FQ: `1.2.a --> 1.2.b`

## 2020-04-10 edition.

- combinator-FQ: calculating of expected genome length is embedded once again -- not considering contigs multiplicity;
- fasta-GC-content is added;

### Version changes:

1. combinator-FQ: `1.1.c --> 1.2.a`

## 2020-03-18 edition.

- combinator-FQ: calculating of expected genome length is disabled because it's impossible to handle high-copy replicons properly;

### Version changes:

1. combinator-FQ: `1.1.b --> 1.1.c`

## 2020-03-09 edition.

- combinator-FQ: calculation of expected genome length fixed;

### Version changes:

1. combinator-FQ: `1.1.a --> 1.1.b`

## 2020-03-07 edition.

- combinator-FQ: output file separated into three: `adjacent_contigs`. `full_matching_log` and `summary`;
- combinator-FQ: A5-coverage summary bug fixed;

### Version changes:

1. combinator-FQ: `1.0.a --> 1.1.a`

## 2020-03-06 edition.

- combinator-FQ version `1.0.a` added;
