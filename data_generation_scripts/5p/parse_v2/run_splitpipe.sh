#!/bin/bash
# Load necessary modules
module load splitpipe
module load Anaconda3
module load STAR
module load SAMtools
module load FastQC

# Define variables for file paths and parameters
genome_dir="/fh/fast/_IRC/FHIL/pub/ref/parse/genomev1.3.1/hg38_V1.3.1"
# fq1="/fh/fast/_IRC/FHIL/grp/BM_Paper_Final_Data/B2-Parse_25k/Subsampled25K_Fastqs/Sublib4/B2-PA-WT-4_S4_R1_292M.fastq.gz"
# fq2="/fh/fast/_IRC/FHIL/grp/BM_Paper_Final_Data/B2-Parse_25k/Subsampled25K_Fastqs/Sublib4/B2-PA-WT-4_S4_R2_292M.fastq.gz"
# output_dir="/fh/fast/_IRC/FHIL/grp/BM_Paper_Final_Data/B2-Parse_25k/OutPut_25K_V.1.1/Sublib4"

SLURM_OPTS="--job-name=satcurveparse \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=8 \
             --mem-per-cpu=10G \
             --time=24:00:00"

samples=(B1-PA-1 B1-PA-2 B1-PA-3 B1-PA-4 B1-PA-5 B1-PA-6 B1-PA-7 B1-PA-8)
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)

# declare -A sample_cell_counts 
# sample_cell_counts=(["F1"]=7098 ["F2"]=9723)
# total_reads=(["F1A"]=411350000 ["F1B"]=412625000 ["F5A"]=328950000 ["F5B"]=336875000)
for target in ${targets[@]}; do
    outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/parse_v2/sat_curves/${target}"
    mkdir $outdir;
    for sample in ${samples[@]}; do
        outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/parse_v2/sat_curves/${target}/${sample}"
        mkdir $outdir;
        cd $outdir;
        # reads=$(echo "scale=5; $target * ${sample_cell_counts[$sample]}" | bc);
        # echo ${sample_cell_counts[${sample}]};
        fq1="/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm01_tcr_5p/parse_v2/sat_curves/${target}/${sample}/${sample}_downsampled_${target}_R1.fastq.gz"
        fq2="/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm01_tcr_5p/parse_v2/sat_curves/${target}/${sample}/${sample}_downsampled_${target}_R2.fastq.gz"
        # sed -e "s/{{SAMPLE}}/${sample}/g" -e "s/{{TARGET}}/$target/g" -e "s/{{CELLS}}/${sample_cell_counts[$sample]}/g" /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm01_tcr_5p/nextgem_5p/sat_curves/config_manual_downsample.csv > $outdir/config.csv
        sbatch $SLURM_OPTS \
        --wrap "split-pipe --mode all \
           --chemistry v2 \
           --genome_dir "$genome_dir" \
           --fq1 "$fq1" \
           --fq2 "$fq2" \
           --output_dir "$outdir" \
           --sample B2-PA-F1A A1-A12 \
           --sample B2-PA-F1B B1-B12 \
           --sample B2-PA-F5A C1-C12 \
           --sample B2-PA-F5B D1-D12 \
           --no_save_anndata --nthreads 8";
    done
done