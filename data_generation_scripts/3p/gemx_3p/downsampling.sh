#!/bin/bash

# Calculate the number of reads per file pair
reads_per_file_pair=$((823200000 / 4))

module load seqtk 
# Loop over each lane and subsample the reads
for lane in {L001,L002,L003,L004}; do
    # Subsample the R1 reads
    seqtk sample -s100 GEM-X_B2-3X-F5A_S7_${lane}_R1_001.fastq.gz $reads_per_file_pair | gzip > subsampled_GEM-X_B2-3X-F5A_S7_${lane}_R1_001.fastq.gz
    
    # Subsample the R2 reads
    seqtk sample -s100 GEM-X_B2-3X-F5A_S7_${lane}_R2_001.fastq.gz $reads_per_file_pair | gzip > subsampled_GEM-X_B2-3X-F5A_S7_${lane}_R2_001.fastq.gz
done