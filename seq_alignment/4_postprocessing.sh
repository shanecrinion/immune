#!/bin/bash

# Directory containing SAM files
sam_dir="/mnt/data/seq_alignment"

# Directory to store processed BAM files
bam_dir="/mnt/data/alignment_postproc/bam_files"

# Iterate over SAM files
for i in {001..200}; do
    sam_file="${sam_dir}/PSR${i}_aligned_reads.sam"
    bam_file="${bam_dir}/PSR${i}_aligned_reads.bam"
    sorted_bam_file="${bam_dir}/PSR${i}_aligned_reads_sorted.bam"
    indexed_bam_file="${bam_dir}/PSR${i}_aligned_reads_sorted.bam.bai"
    dedup_bam_file="${bam_dir}/PSR${i}_aligned_reads_sorted_dedup.bam"
    dedup_metrics="${bam_dir}/PSR${i}_aligned_reads_dedup_metrics.txt"

    # Sort SAM file
    samtools sort -o "$sorted_bam_file" "$sam_file"

    # Convert to BAM format
    samtools view -b -o "$bam_file" "$sorted_bam_file"

    # Index the BAM file
    samtools index "$bam_file" "$indexed_bam_file"

    # Remove PCR duplicates
    samtools rmdup -s "$sorted_bam_file" "$dedup_bam_file" > "$dedup_metrics"

    # Optional: Generate QC metrics
    samtools flagstat "$dedup_bam_file" > "${bam_dir}/PSR${i}_aligned_reads_flagstat.txt"
done

echo "Processing complete."

