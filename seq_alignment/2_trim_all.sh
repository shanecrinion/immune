#!/bin/bash

# Loop through all files with the prefix "PSR" and the suffix "_R1_001.fastq.gz"
for file in PSR*_R1_001.fastq.gz; do
        # Extract the base filename without the extension and "_R1_001"
	filename=$(basename "$file" _R1_001.fastq.gz)
        # Perform adaptor trimming using Trimmomatic for the corresponding R1 and R2 files
        java -jar /home/administrator/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
            "${filename}_R1_001.fastq.gz" "${filename}_R2_001.fastq.gz" \
            "${filename}_R1_001_TRIM_paired.fastq.gz" "${filename}_R1_001_TRIM_unpaired.fastq.gz" \
            "${filename}_R2_001_TRIM_paired.fastq.gz" "${filename}_R2_001_TRIM_unpaired.fastq.gz" \
            ILLUMINACLIP:/home/administrator/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done


