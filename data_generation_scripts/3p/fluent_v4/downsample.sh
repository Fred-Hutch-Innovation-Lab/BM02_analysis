#!/bin/bash

# Job Options- must be before *any* executable lines

#SBATCH --job-name="downsample"
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=36       # cpu-cores per task (>1 if multi-threaded tasks) ## Gizmo nodes have 36 cores
#SBATCH --mem-per-cpu=8G        # memory per cpu-core
#SBATCH --time=12:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=dgratz@fredhutch.org

module load seqtk/1.3-GCC-10.2.0
module load pigz/2.8-GCCcore-13.2.0

fastq_dir="/fh/fast/_IRC/FHIL/grp/BM02_FluentV4_2024/fastqs/full"
outdir="/fh/fast/_IRC/FHIL/grp/BM02_FluentV4_2024/fastqs/downsampled"

declare -A target_reads

# target_reads["F1A"]=449738139 # Downsampling not needed
# target_reads["F1B"]=415600000
target_reads["F5A"]=419575000
# target_reads["F5B"]=380850000

for sample in "${!target_reads[@]}"; do
    # echo $sample ${target_reads[$sample]};
    for filetype in R1 R2; do
        seqtk sample -s100 ${fastq_dir}/*${sample}*${filetype}*.fastq.gz ${target_reads[$sample]} | pigz -p 36 > ${outdir}/B2-FBv4_${sample}_downsampled_${filetype}.fastq.gz
    done
done
