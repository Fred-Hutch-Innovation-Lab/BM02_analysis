#!/bin/bash

# Load seqtk module
module load seqtk

# Change directory
cd /fh/fast/_IRC/FHIL/grp/BM01_5Data10X/FH_NextGEM_5p_fastqs/F1/GEX

# Calculate the number of reads per file pair
reads_per_file_pair=$((179900000 / 4))

# Loop through each lane and perform subsampling
for lane in L001 L002 L003 L004; do
    # Subsample R1 files and compress the output
    seqtk sample -s100 NextGEM_B1_5P_F1_gex_S1_${lane}_R1_001.fastq.gz $reads_per_file_pair | gzip > subsampled_NextGEM_B1_5P_F1_gex_S1_${lane}_R1_001.fastq.gz
    
    # Subsample R2 files and compress the output
    seqtk sample -s100 NextGEM_B1_5P_F1_gex_S1_${lane}_R2_001.fastq.gz $reads_per_file_pair | gzip > subsampled_NextGEM_B1_5P_F1_gex_S1_${lane}_R2_001.fastq.gz
done

