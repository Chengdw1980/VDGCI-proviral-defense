#!/bin/bash
# ============================================================================
# VDGCI pipeline — configuration (Anvil-ready, v2, 2026-10)
# Every script sources this file. On Purdue Anvil the defaults below already
# point at the right storage; on another cluster edit the three *_BASE lines.
# ============================================================================
export PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# --- Your Anvil CPU allocation account (run `mybalance` on Anvil to see it,
#     e.g. something like  bio260593-cpu ). REQUIRED for sbatch. ---
export ALLOC="${ALLOC:-bio260593}"   # confirmed with mybalance on 2026-10-23

# --- Which genome list to run (one NCBI accession per line, GCF_/GCA_) ---
#   config/accessions_test10.txt   10 control genomes  (run this FIRST)
#   config/accessions_HQ_5054.txt  primary cohort      (GTDB r232, 1 per family, HQ)
#   config/accessions_full_6048.txt sensitivity cohort  (all named families)
#   config/accessions_48.txt       the published 48-genome pilot
export ACCESSION_LIST="${ACCESSION_LIST:-${PIPELINE_ROOT}/config/accessions_test10.txt}"
export RUN_NAME="${RUN_NAME:-$(basename "$ACCESSION_LIST" .txt | sed 's/^accessions_//')}"

# --- Storage. Anvil: $SCRATCH = fast, 100 TB, PURGED after 30 days unused;
#     $PROJECT = 5 TB, kept while the allocation is active (snapshotted). ---
SCRATCH_BASE="${SCRATCH:-${PIPELINE_ROOT}/scratch}"
PROJECT_BASE="${PROJECT:-${PIPELINE_ROOT}/project}"
export WORK="${WORK_OVERRIDE:-${SCRATCH_BASE}/vdgci/${RUN_NAME}}"      # big, temporary
export KEEP="${KEEP_OVERRIDE:-${PROJECT_BASE}/vdgci/${RUN_NAME}}"      # results kept safe
export DB_DIR="${DB_DIR:-${PROJECT_BASE}/vdgci/db}"                   # databases: build once, keep
export PROJECT_BASE
export CONDA_ENV="${CONDA_ENV:-${PROJECT_BASE}/vdgci/env}"
export TOOLS_BIN="${TOOLS_BIN:-${PROJECT_BASE}/vdgci/bin}"   # v3.1: pinned diamond + defense-finder wrappers         # conda env lives in project space (home is only 25 GB)

export GENOMES_DIR="${WORK}/genomes"     # <ACC>.fna
export PROT_DIR="${WORK}/prot_by_genome" # <ACC>.faa  (NCBI proteins, if annotated)
export GBFF_DIR="${WORK}/gbff"           # <ACC>.gbff (only for NCBI-annotated genomes)
export GENOMAD_OUT="${WORK}/genomad_out"
export DF_OUT="${WORK}/defensefinder_out"
export PHISPY_OUT="${WORK}/phispy_out"
export CHECKV_OUT="${WORK}/checkv_out"
export PHYLO_OUT="${WORK}/phylo"
export LOGS="${WORK}/logs"
export RESULTS="${KEEP}/results"         # final tables go straight to safe storage

export GENOMAD_DB="${DB_DIR}/genomad_db"
export CHECKV_DB="${DB_DIR}/checkv-db-v1.5"
export DF_MODELS="${DB_DIR}/df_models"
export GTDB_TREE="${DB_DIR}/gtdb/bac120_r232.tree"

# --- Per-genome compute. Anvil regular node = 128 cores, ~257 GB (~2 GB/core).
#     8 threads x 16 GB per genome = 16 genomes packed per 128-core node. ---
export THREADS="${THREADS:-8}"
export MEM_PER_GENOME_GB="${MEM_PER_GENOME_GB:-16}"
export GENOMAD_SPLITS="${GENOMAD_SPLITS:-8}"   # keeps geNomad's marker search inside 16 GB
export GENOMES_PER_NODE="${GENOMES_PER_NODE:-16}"

export GENOMAD_MINSCORE="0.7"            # geNomad default; recorded for provenance
export CONSENSUS_OVERLAP="0.5"           # geNomad<->PhiSpy reciprocal-overlap threshold

mkdir -p "$WORK" "$KEEP" "$DB_DIR" "$GENOMES_DIR" "$PROT_DIR" "$GBFF_DIR" "$GENOMAD_OUT" \
         "$DF_OUT" "$PHISPY_OUT" "$CHECKV_OUT" "$PHYLO_OUT" "$LOGS" "$RESULTS" 2>/dev/null
export N_GENOMES=$(grep -cve '^\s*$' "$ACCESSION_LIST" 2>/dev/null || echo 0)
