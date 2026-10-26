#!/bin/bash
# Run the whole pipeline serially on one machine (small sets only, e.g. the 10-genome test).
set -euo pipefail
cd "$(dirname "$0")/.."; source config/config.sh
bash scripts/01_setup_databases.sh
bash scripts/02_fetch_genomes.sh
while read -r A; do [ -z "$A" ] || bash scripts/03_run_per_genome.sh "$A"; done < "$ACCESSION_LIST"
bash scripts/04_phylogeny.sh || echo "(16S tree step skipped)"
python3 scripts/05_analyze.py
python3 scripts/07_cost_report.py "$N_GENOMES"
