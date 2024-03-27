#!/usr/bin/awk -f
#
# Script for summarizing variation ant single specified position in SAM/BAM file.
# Tested on GNU Awk 5.0.1.
#
# Version 0.2.a
# 2024-03-27
#
# == Usage ==
# 1
# Basic usage:
#    samtools view <BAM/SAM> | awk -f sum_up_snv.awk pos=<POS>
# Example:
#    samtools view my_mapping.bam | awk -f sum_up_snv.awk pos=1067
# Verbose: for each read print a tab-separated table having three columns: read_name, base, base_quality
#    samtools view my_mapping.bam | awk -f sum_up_snv.awk pos=1067 verbose=1
# 2
# Faster -- specify reference sequence name:
#    samtools view <BAM/SAM> <RNAME> | awk -f sum_up_snv.awk pos=<POS>
# Example:
#    samtools view my_mapping.bam chr1 | awk -f sum_up_snv.awk pos=1067
# 3
# Even faster -- specify single-base region:
#    samtools view <BAM/SAM> <RNAME>:<POS>-<POS> | awk -f sum_up_snv.awk pos=<POS>
# Example:
#    samtools view my_mapping.bam chr1:1067-1067 | awk -f sum_up_snv.awk pos=1067


BEGIN{
    cigar_operations="MIDHS";
    a_count = 0;
    t_count = 0;
    g_count = 0;
    c_count = 0;
    del_count = 0;
    for (ascii_code = 0; ascii_code < 256; ascii_code++) {
        ord[sprintf("%c", ascii_code)] = ascii_code;
    }
}

{
    # $4 -- alignment start pos; $10 -- read sequence.
    # If read does not map to `pos` -- omit it immidiately.
    if ($4 > pos || $4 + length($10) <= pos) {
        next; # go to next line in bam/sam file
    }

    start_aln = $4;
    seq = toupper($10);

    cigar = $6;
    len_cigar = length(cigar);

    # Position "in reference axis"
    ref_pos = start_aln - 1;
    # Position "in read axis"
    read_pos = 0;

    # Current CIGAR operation, e.g. '163M'
    curr_cirar = "";
    # Length of `curr_cirar`, e.g. '163'
    curr_operation_len = "";

    # Go through CIGAR string in order to take into account all indels and hard/soft clips
    for (i = 1; i <= len_cigar; i++) {

        # Current character of CIGAR string
        chr = substr(cigar, i, 1);
        # Update CIGAR
        curr_cirar = sprintf("%s%s", curr_cirar, chr);

        # If `chr` is a letter -- it's a CIGAR operation
        if ( index(cigar_operations, chr) != 0 ) {

            curr_operation_len = substr(curr_cirar, 1, length(curr_cirar)-1);

            # Go further "in reference axis"
            if (chr == "M" || chr == "D") {
                ref_pos += curr_operation_len;
            }
            # Go further "in read axis"
            if (chr != "D" && chr != "H") {
                read_pos += curr_operation_len;
            }

            # If we reach our `pos`
            if (ref_pos >= pos) {

                # Not a deletion (not "*") but a base
                if (chr != "D") {

                    base_pos = read_pos - (ref_pos-pos);
                    base = substr(seq, base_pos, 1);

                    if (verbose) {
                        qual_encoded = substr($11, base_pos, 1);
                        qual = ord[qual_encoded] - 33;
                        printf("%s\t%s\t%d\n", $1, base, qual);
                    }

                    if (base == "A") {
                        a_count++;
                    } else if (base == "T") {
                        t_count++;
                    } else if (base == "G") {
                        g_count++;
                    } else if (base == "C") {
                        c_count++;
                    }
                } else {
                    # Deletion -- "*"
                    del_count++;
                    if (verbose) {
                        printf("%s\t*\t-\n", $1)
                    }
                }
                break;
            }
            # Clear CIGAR
            curr_cirar = "";
            curr_operation_len = "";
        }
    }
}

END{
    cov = a_count + t_count + g_count + c_count + del_count;
    printf("Position %d:\n", pos)
    printf("Coverage: %d\n", cov)

    cov_ratio = 100 / cov; 
    printf("A:%5d | %5.1f%%\n",   a_count,   a_count * cov_ratio);
    printf("T:%5d | %5.1f%%\n",   t_count,   t_count * cov_ratio);
    printf("G:%5d | %5.1f%%\n",   g_count,   g_count * cov_ratio);
    printf("C:%5d | %5.1f%%\n",   c_count,   c_count * cov_ratio);
    printf("*:%5d | %5.1f%%\n", del_count, del_count * cov_ratio);
}
