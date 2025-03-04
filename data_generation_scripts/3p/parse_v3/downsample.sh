#!/bin/bash

SLURM_OPTS="--job-name=downsampleparse \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=4 \
             --mem-per-cpu=16G \
             --time=40:00:00 \
             --mail-type=begin \
             --mail-type=end \
             --mail-type=fail \
             --mail-user=dgratz@fredhutch.org"

fastq_dir=/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/parse_v3/downsampled/
outdir=/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/parse_v3/sat_curves/
# for sublib_dir in $fastq_dir/Sublib*; do
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000)

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

    target_total_reads=$(echo "scale=5; $target_reads / 25000" | bc)
    # local target_total_reads=$((target_reads / 25000))
    echo "Downsampling $sample to $target_total_reads reads"
    for file in R1 R2; do
        # echo $file
        input="${fastq_dir}/${sample}/*${file}*.fastq.gz"
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
    for sample in Sublib1 Sublib2 Sublib3 Sublib4 Sublib5 Sublib6 Sublib7 Sublib8; do
        # echo "Downsampling $sample"
        # mkdir "${outdir}/${target}"
        writedir="${outdir}/${target}/${sample}"
        # mkdir "$writedir"
        # total_reads="${sample_total_reads[$sample]}"
        # cells="${sample_cell_counts[$sample]}"
        downsample "$sample" "$target" "$writedir" #files[@] "$total_reads" "$cells"  & 
    done
    # wait
done
