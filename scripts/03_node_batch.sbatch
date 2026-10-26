#!/bin/bash
# PRODUCTION RUN: one whole 128-core Anvil node processes a CHUNK of genomes,
# 16 at a time (8 threads x 16 GB each), until the chunk is finished.
# This avoids Anvil's per-user job limits (a 5,000-task array would be refused)
# and keeps the node fully busy. Charged 128 SUs per node-hour.
# Do not submit by hand — use:  bash scripts/03_submit_chunks.sh
#SBATCH --job-name=vdgci
#SBATCH --partition=wholenode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=48:00:00
#SBATCH --output=logs_slurm/chunk_%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
source config/config.sh
module load anaconda
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"
CHUNK="$1"
echo "Chunk $CHUNK: $(wc -l < "$CHUNK") genomes, ${GENOMES_PER_NODE} in parallel, start $(date -Is)"
xargs -a "$CHUNK" -P "$GENOMES_PER_NODE" -I{} bash scripts/03_run_per_genome.sh {}
echo "Chunk $CHUNK finished $(date -Is)"
