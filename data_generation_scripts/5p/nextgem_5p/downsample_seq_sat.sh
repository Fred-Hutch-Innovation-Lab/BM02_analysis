#!/bin/bash

# Job Options- must be before *any* executable lines

module load seqtk/1.3-GCC-10.2.0
module load pigz/2.8-GCCcore-13.2.0

sampledir='/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm01_tcr_5p/nextgem_5p/downsampled/Subsampled25k_WT_Fastqs'
outdir='/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm01_tcr_5p/nextgem_5p/sat_curves/'
targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000) #x

declare -A sample_cell_counts 
sample_cell_counts=(["F1"]=7098 ["F2"]=9723)
sample_total_reads=(["F1"]=179900000 ["F2"]=245175000)

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
    echo $target 
    echo $cells
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
        # mkdir "${outdir}/${sample}"
        writedir="${outdir}/${sample}/${target}"
        # mkdir "$writedir"
        total_reads="${sample_total_reads[$sample]}"
        cells="${sample_cell_counts[$sample]}"
        downsample "$sample" "$total_reads" "$cells" "$target" "$writedir" #files[@] &
    done
    # wait
done
