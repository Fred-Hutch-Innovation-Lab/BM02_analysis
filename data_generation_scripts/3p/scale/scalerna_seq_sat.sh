#!/bin/bash

# Job Options- must be before *any* executable lines

#SBATCH --job-name="satcurvescale"
#SBATCH --output=sbatch_scale.out
#SBATCH --error=sbatch_scale.err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16        # cpu-cores per task (>1 if multi-threaded tasks) ## Gizmo nodes have 36 cores
#SBATCH --mem-per-cpu=12G        # memory per cpu-core
#SBATCH --time=48:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=dgratz@fredhutch.org

module load Java
module load Singularity/3.5.3

cd /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/scale/sat_curves

targets=(1000 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000) # ) #x
datadir="/fh/fast/_IRC/FHIL/grp/BM_paper/fastqs/bm02_wt_3p/scale/sat_curves"
outdir="/fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/scale/sat_curves"
samples=(F1A F1B F5A F5B)

SLURM_OPTS="--job-name=satcurvescale \
             --output=sbatch.out \
             --error=sbatch.err \
             --nodes=1 \
             --ntasks=1 \
             --cpus-per-task=8 \
             --mem-per-cpu=8G \
             --time=24:00:00 \
             --mail-type=begin \
             --mail-type=end \
             --mail-type=fail \
             --mail-user=dgratz@fredhutch.org"

for target in "${targets[@]}"; do
    mkdir ${outdir}/${target}
    for sample in "${samples[@]}"; do
        for file in ${datadir}/${sample}/${target}/*.fastq.gz; do
            ln -s -T $file ${outdir}/${target}/seqsat_$(basename ${file})
            # ln -s -T $file ${outdir}/${target}/seq.sat_${file}
        done
    done
        # ln -s ${datadir}/${sample}/*.fastq.gz ${outdir}/${target}/
    
    sbatch $SLURM_OPTS \
        --wrap "
        /fh/fast/_IRC/FHIL/grp/bioinfo_tools/nextflow/nextflow run \
        /fh/fast/_IRC/FHIL/grp/bioinfo_tools/ScaleRNA/ScaleRna/ \
        -profile singularity \
        -params-file /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/scale/sat_curves/runParams.yml \
        --fastqDir ${outdir}/${target} \
        --outDir ${outdir}/${target}/scaleRNA.out
        ";
done
