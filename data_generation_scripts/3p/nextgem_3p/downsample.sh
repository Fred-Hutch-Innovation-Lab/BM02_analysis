#!/bin/bash

# Load the seqtk module
module load seqtk

# Calculate the number of reads per file pair
reads_per_file_pair=$((411350000 / 4))

# Loop over each lane and subsample the reads
for lane in {L001,L002,L003,L004}; do
    # Subsample the R1 reads
    seqtk sample -s100 /NextGEM_B2-3P-F1A_S1_${lane}_R1_001.fastq.gz $reads_per_file_pair | gzip > subsampled_NextGEM_B2-3P-F1A_S1_${lane}_R1_001.fastq.gz
    
    # Subsample the R2 reads
    seqtk sample -s100 NextGEM_B2-3P-F1A_S1_${lane}_R2_001.fastq.gz $reads_per_file_pair | gzip > subsampled_NextGEM_B2-3P-F1A_S1_${lane}_R2_001.fastq.gz
done
