#!/bin/bash
# SMALL TEST RUN (the 10 control genomes) on Anvil's 'shared' partition.
# One genome per array task, 8 cores + 16 GB each -> charged 8 SUs per hour.
# Submit from the pipeline folder:
#   sbatch -A <your-alloc> scripts/03_test_shared.sbatch
#SBATCH --job-name=vdgci-test
#SBATCH --partition=shared
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs_slurm/test_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export ACCESSION_LIST="${ACCESSION_LIST:-$PWD/config/accessions_test10.txt}"
source config/config.sh
module load anaconda
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSION_LIST")
[ -n "$ACC" ] && bash scripts/03_run_per_genome.sh "$ACC"
