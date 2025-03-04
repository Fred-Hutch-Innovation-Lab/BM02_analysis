#!/bin/bash
# Load necessary modules
module load splitpipe
module load STAR
module load SAMtools
module load FastQC

# Define variables for file paths and parameters
genome_dir="/fh/fast/_IRC/FHIL/pub/ref/parse/genomev1.1.3/hg38_V1.1.3"
fq1="/fh/fast/_IRC/FHIL/grp/BM_Paper_Final_Data/B2-Parse_25k/Subsampled25K_Fastqs/Sublib4/B2-PA-WT-4_S4_R1_292M.fastq.gz"
fq2="/fh/fast/_IRC/FHIL/grp/BM_Paper_Final_Data/B2-Parse_25k/Subsampled25K_Fastqs/Sublib4/B2-PA-WT-4_S4_R2_292M.fastq.gz"
output_dir="/fh/fast/_IRC/FHIL/grp/BM_Paper_Final_Data/B2-Parse_25k/OutPut_25K_V.1.1/Sublib4"

# Run split-pipe with the specified parameters
split-pipe --mode all \
           --chemistry v3 \
           --genome_dir "$genome_dir" \
           --fq1 "$fq1" \
           --fq2 "$fq2" \
           --output_dir "$output_dir" \
           --sample B2-PA-F1A A1-A12 \
           --sample B2-PA-F1B B1-B12 \
           --sample B2-PA-F5A C1-C12 \
           --sample B2-PA-F5B D1-D12
