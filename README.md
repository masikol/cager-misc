# CAGER-misc

Miscellaneous bioinformatic tools from [Laboratory](https://mbio.bas-net.by/cager/en/) “The **C**enter of **A**nalytical and **G**enetic **E**ngineering **R**esearch”.

List of tools with links to manuals:
- [sum-up-snv](https://www.github.com/masikol/cager-misc/wiki/sum-up-snv): a script for counting coverage and single nucleotide variants at a single specified position in SAM/BAM file;
- [cigar_maplen](https://www.github.com/masikol/cager-misc/wiki/cigar_maplen): a script for manual inspection of structural variants, mainly for long reads. It prints which portion of a read is actually mapped, and which portions are clipped from each side;
- [samtools_setop](https://www.github.com/masikol/cager-misc/wiki/samtools_setop): the script performs set operations (intersection, union, difference) on read IDs that are mapped to arbitrary reference positions;
- [gla-glar](https://www.github.com/masikol/cager-misc/wiki/gla-glar): “GLAde GLARer”: the script finds long intergenic and inter-CDS regions in GenBank files;
- AAI-taxon.sh(Help page is yet to be done): the script calculates Average Amino acid Identities between a query genome and all type strain genomes of a given taxon;
- [pub](https://www.github.com/masikol/cager-misc/wiki/pub): a script for automatic selection of sequencing barcodes;
- [dedupl-fastq](https://www.github.com/masikol/cager-misc/wiki/dedupl-fastq): the script is designed for deduplication of fastq files;
- [mean-qual](https://www.github.com/masikol/cager-misc/wiki/mean-qual): the script calculates mean quality of reads in `fastq` file(s);
- [most-freq-subseq](https://www.github.com/masikol/cager-misc/wiki/most-freq-subseq): the script finds N most frequently occuring subsequences of given length for each sequence in fasta file;
- [NOS](https://www.github.com/masikol/cager-misc/wiki/NOS): the sript counts non-overalapping occurences of query sequence (and it's reverse complement "comrade") in `fasta` file(s);
- [fasta-GC-content](https://www.github.com/masikol/cager-misc/wiki/fasta-GC-content): the script calculates GC-content of each sequence in `fasta` file(s);
- [fastq2fasta](https://www.github.com/masikol/cager-misc/wiki/fastq2fasta): the script converts `fastq` files to `fasta` format;
- [fastq-read-count](https://www.github.com/masikol/cager-misc/wiki/fastq-read-count): the script counts amount of reads in `fastq` file(s);
- [find-seq](https://www.github.com/masikol/cager-misc/wiki/find-seq): the script finds fasta record(s) in `fasta` file by given sequence header;
- [dna-summary](https://www.github.com/masikol/cager-misc/wiki/dna-summary): the script collects basic information from `.dna` [SPAdes](http://cab.spbu.ru/software/spades/) contigs in `contigs/` directory;
- [packer-dna-to-fasta](https://www.github.com/masikol/cager-misc/wiki/packer-dna-to-fasta): the script packs `.dna` [SPAdes](http://cab.spbu.ru/software/spades/) contigs in 'contigs' directory to single multi-fasta file;
- [seqator](https://www.github.com/masikol/cager-misc/wiki/seqator): the script moves `.dna` [SPAdes](http://cab.spbu.ru/software/spades/) contigs with coverage less than specified one from `contigs/` directory to directory `cov_below_x/`;
- combinator-FQ: genome assembly facilitation. This script is now moved to the separate repository: [https://github.com/masikol/combinator-FQ](https://github.com/masikol/combinator-FQ);
- kromsatel: a tool for splitting chimeric nanopore amplicon reads. This script is now moved to the separate repository: [https://github.com/masikol/kromsatel](https://github.com/masikol/kromsatel);
