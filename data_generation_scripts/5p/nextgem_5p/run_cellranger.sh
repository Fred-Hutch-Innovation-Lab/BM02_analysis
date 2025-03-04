#!/bin/bash
#SBATCH --job-name=cellranger_run     # Job name
#SBATCH --output=cellranger_run.out   # Output file for the job
#SBATCH --error=cellranger_run.err    # Error file for the job
#SBATCH --time=24:00:00               # Maximum time the job can run (e.g., 24 hours)
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=1                    # Number of tasks (usually 1 for CellRanger)
#SBATCH --cpus-per-task=16            # Number of CPUs per task (adjust based on your needs)
#SBATCH --mem=64G                     # Memory needed for the job

# Load the CellRanger module
module load CellRanger

# Change to the appropriate directory
cd /fh/fast/_IRC/FHIL/grp/BM_Paper_Final_Data/B1-NextGEM5P/

# Run CellRanger multi with the specified configuration
cellranger multi --id=NextGEM_F1 --csv=./config.csv

