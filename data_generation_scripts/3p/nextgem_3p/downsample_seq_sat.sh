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

sampledir='/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/nextgem_3p/downsampled/'
outdir='/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/nextgem_3p/sat_curves/'
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000) #x

declare -A sample_cell_counts 
sample_cell_counts=(["F1A"]=16338 ["F1B"]=16403 ["F5A"]=13078 ["F5B"]=13413)
sample_total_reads=(["F1A"]=411350000 ["F1B"]=412625000 ["F5A"]=328950000 ["F5B"]=336875000)

SLURM_OPTS="--job-name=downsample_nextgem \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=2 \
             --mem-per-cpu=16G \
             --time=1:00:00 
             "

downsample() {
    local sample=$1
    local total_reads=$2
    local cells=$3
    local target_reads=$4
    local writedir=$5

    local target_total_reads=$((target_reads * cells / 4))
    echo "Downsampling $sample to $target_total_reads reads"
    for file in R1 R2; do #"${files[@]}"; do
        echo $file
        for lane in L001 L002 L003 L004; do
            input="${sampledir}/${sample}/*${sample}*${lane}*${file}*.fastq.gz"
            sbatch $SLURM_OPTS \
                    --wrap "seqtk sample \
                                -s100 $input $target_total_reads | \
                                pigz > ${writedir}/${sample}_downsampled_${target_reads}_S1_${lane}_${file}_001.fastq.gz";
        done
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
