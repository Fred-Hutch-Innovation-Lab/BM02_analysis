#!/bin/bash

SLURM_OPTS="--job-name=pipseeker \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=16 \
             --mem-per-cpu=10G \
             --time=24:00:00 \
             --mail-type=begin \
             --mail-type=end \
             --mail-type=fail \
             --mail-user=dgratz@fredhutch.org"

samples=(B2-FBv5_F1A B2-FBv5_F1B B2-FBv5_F5A B2-FBv5_F5B)

for sample in ${samples[@]}; do
    outdir="/fh/fast/_IRC/FHIL/grp/BM02_FluentV5/downsampled_runs/${sample}"
    cd $outdir;
    sbatch $SLURM_OPTS \
    --wrap "/fh/fast/_IRC/FHIL/grp/bioinfo_tools/pipseeker/pipseeker-v3.3.0-linux/pipseeker full \
            --fastq /fh/fast/_IRC/FHIL/grp/BM02_FluentV5/fastqs/downsampled/${sample} \
            --star-index-path /shared/biodata/ngs/Reference/PIPseeker/pipseeker-gex-reference-GRCh38-2022.04 \
            --chemistry V \
            --annotation /shared/biodata/ngs/Reference/PIPseeker/human-pbmc-references/references/human-pbmc-v4-detailed.csv \
            --output-path ${outdir}";
done