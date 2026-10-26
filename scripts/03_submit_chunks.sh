#!/bin/bash
# Split the not-yet-finished genomes into chunks and submit one whole-node job
# per chunk. Re-running it later only resubmits genomes that are not done.
#   Usage:  bash scripts/03_submit_chunks.sh [genomes_per_chunk]   (default 320)
# 320 genomes / 16 at a time = 20 "waves"; at <=2 h per wave that fits the 48 h limit.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
[ "$ALLOC" = "CHANGE_ME-cpu" ] && { echo "Set ALLOC first (see config/config.sh; run mybalance)"; exit 1; }
PER="${1:-320}"; mkdir -p logs_slurm "${KEEP}/chunks"
TODO="${KEEP}/chunks/todo_$(date +%Y%m%d_%H%M).txt"; : > "$TODO"
while read -r A; do [ -z "$A" ] || [ -s "${KEEP}/done/${A}" ] || echo "$A" >> "$TODO"; done < "$ACCESSION_LIST"
N=$(wc -l < "$TODO"); echo "${N} genomes to run (of ${N_GENOMES}); ${PER} per node-job"
[ "$N" -gt 0 ] || exit 0
split -l "$PER" -d -a 3 "$TODO" "${TODO%.txt}_chunk_"
for C in "${TODO%.txt}"_chunk_*; do
  sbatch -A "$ALLOC" scripts/03_node_batch.sbatch "$C"
done
echo "Submitted. Watch with:  squeue -u \$USER    and   ls ${KEEP}/done | wc -l"
