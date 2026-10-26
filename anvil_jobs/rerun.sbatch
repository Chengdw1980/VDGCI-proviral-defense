#!/bin/bash
#SBATCH --job-name=vdgci-dedup
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/home/x-dcheng3/rerun_dedup_%j.out
set -uo pipefail
cd /anvil/projects/x-bio260593/VDGCI_cluster_pipeline
unset RUN_NAME KEEP WORK RESULTS
export ACCESSION_LIST=$PWD/config/accessions_HQ_5054.txt
source config/config.sh
module load anaconda; eval "$(conda shell.bash hook)"; conda activate "$CONDA_ENV"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
echo "== 05 core $(date -Is)"; python3 scripts/05_analyze.py
echo "== 06 extended $(date -Is)"; python3 scripts/06_extended_analysis.py
PY=/anvil/projects/x-bio260593/vdgci/env/bin/python
echo "== defmap $(date -Is)"; $PY ~/defmap.py && $PY ~/defstats.py && $PY ~/defmodel.py && $PY ~/defaudit.py && $PY ~/defrobust.py
echo "== ALL DONE $(date -Is)"
