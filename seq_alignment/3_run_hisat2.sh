#!/bin/bash

# Define the range of file numbers
for i in {001..200}; do
    # Define the file prefix
    prefix="PSR${i}"

    # Run HISAT2 for each pair of files
    hisat2 -x /mnt/data/seq_alignment/grch37_snp_tran/genome_snp_tran -1 "${prefix}_R1_001_TRIM_paired.fastq.gz" -2 "${prefix}_R2_001_TRIM_paired.fastq.gz" -S "${prefix}_aligned_reads.sam"
done

