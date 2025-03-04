#!/bin/bash

# Job Options- must be before *any* executable lines

#SBATCH --job-name="scaleWcells"
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=32        # cpu-cores per task (>1 if multi-threaded tasks) ## Gizmo nodes have 36 cores
#SBATCH --mem-per-cpu=10G        # memory per cpu-core
#SBATCH --time=48:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=dgratz@fredhutch.org

module load Java
module load Singularity/3.5.3

cd /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/scale/downsampled_runs/

/fh/fast/_IRC/FHIL/grp/bioinfo_tools/nextflow/nextflow run \
    /fh/fast/_IRC/FHIL/grp/bioinfo_tools/ScaleRNA/ScaleRna/ \
    -profile singularity \
    -params-file /fh/fast/_IRC/FHIL/grp/BM_paper/processing/bm02_wt_3p/scale/downsampled_runs/runParams_w_cell_calling.yml

