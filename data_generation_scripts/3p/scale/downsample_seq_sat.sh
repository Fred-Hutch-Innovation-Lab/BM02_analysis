#!/bin/bash

# Job Options- must be before *any* executable lines

#SBATCH --job-name="scaleSatCurve"
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=36       # cpu-cores per task (>1 if multi-threaded tasks) ## Gizmo nodes have 36 cores
#SBATCH --mem-per-cpu=8G        # memory per cpu-core
#SBATCH --time=120:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=dgratz@fredhutch.org

module load seqtk/1.3-GCC-10.2.0
module load pigz/2.8-GCCcore-13.2.0

sampledir='/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/scale/downsampled/'
outdir='/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/scale/sat_curves/'
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000) #x

declare -A sample_cell_counts
sample_cell_counts=(["F1A"]=33298 ["F1B"]=30025 ["F5A"]=15782 ["F5B"]=39642)
sample_total_reads=(["F1A"]=832638582 ["F1B"]=745904763 ["F5A"]=397218677 ["F5B"]=991075196)

SLURM_OPTS="--job-name=downsample_scale \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=8 \
             --mem-per-cpu=8G \
             --time=1:00:00 \
             --mail-type=begin \
             --mail-type=end \
             --mail-type=fail \
             --mail-user=dgratz@fredhutch.org"

downsample() {
    local sample=$1
    local total_reads=$2
    local cells=$3
    local target_reads=$4
    local writedir=$5

    local target_total_reads=$((target_reads * cells))
    echo "Downsampling $sample to $target_total_reads reads"
    for file in I1 I2 R1 R2; do #"${files[@]}"; do
        echo $file
        input="${sampledir}/*${sample}*${file}*.fastq.gz"
        sbatch $SLURM_OPTS \
                --wrap "seqtk sample \
                            -s100 $input $target_total_reads | \
                            pigz > ${writedir}/${sample}_downsampled_${target_reads}_${file}.fastq.gz";
    done
}

# Iterate over samples
for target in "${targets[@]}"; do
    printf "\n-------------\nDownsampling to ${target} reads per cell"
    for sample in "${!sample_cell_counts[@]}"; do
        echo "Downsampling $sample"
        mkdir "${outdir}/${sample}"
        writedir="${outdir}/${sample}/${target}"
        mkdir "$writedir"
        total_reads="${sample_total_reads[$sample]}"
        cells="${sample_cell_counts[$sample]}"
        downsample "$sample" "$total_reads" "$cells" "$target" "$writedir" #files[@] &
    done
    # wait
done
