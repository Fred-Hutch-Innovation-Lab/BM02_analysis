#!/bin/bash

SLURM_OPTS="--job-name=downsampleparse \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=4 \
             --mem-per-cpu=16G \
             --time=12:00:00
             "

fastq_dir=/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm01_tcr_5p/parse_v2/full/
outdir=/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm01_tcr_5p/parse_v2/sat_curves/
# for sublib_dir in $fastq_dir/Sublib*; do
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500)

# declare -A sample_cell_counts
# sample_cell_counts=(["F1A"]=33298 ["F1B"]=30025 ["F5A"]=15782 ["F5B"]=39642)
# sample_total_reads=(["F1A"]=832638582 ["F1B"]=745904763 ["F5A"]=397218677 ["F5B"]=991075196)

module load seqtk/1.3-GCC-10.2.0
module load pigz/2.8-GCCcore-13.2.0

downsample() {
    local sample=$1
    # local total_reads=$2
    # local cells=$3
    local target_reads=$2
    local writedir=$3

    target_total_reads=$(echo "scale=5; $target_reads * 11780" | bc) #11780 = sum of cells per sample / 8 sublibraries
    echo "Downsampling $sample to $target_total_reads reads"
    for file in R1 R2; do
        # echo $file
        input="${fastq_dir}/${sample}*${file}*.fastq.gz"
        sbatch $SLURM_OPTS \
            --wrap "#module load seqtk/1.3-GCC-10.2.0
                    #module load pigz/2.8-GCCcore-13.2.0
                    seqtk sample -s100 $input $target_total_reads | \
                        pigz > ${writedir}/${sample}_downsampled_${target_reads}_${file}.fastq.gz";
    done
}

# Iterate over samples
for target in "${targets[@]}"; do
    printf "\n-------------\nDownsampling to ${target} reads per cell"
    for sample in B1-PA-1 B1-PA-2 B1-PA-3 B1-PA-4 B1-PA-5 B1-PA-6 B1-PA-7 B1-PA-8; do
        # mkdir "${outdir}/${target}"
        writedir="${outdir}/${target}/${sample}"
        # mkdir "$writedir"
        downsample "$sample" "$target" "$writedir" #files[@] "$total_reads" "$cells"  & 
    done
    # wait
done
