#!/bin/bash
# Download reference databases ONCE into $DB_DIR (project storage, so Anvil's
# scratch purge cannot delete them). Run on an Anvil LOGIN node (needs internet).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
export PATH="${TOOLS_BIN}:${PATH}"   # pinned DIAMOND 2.0.15 + DefenseFinder 3.0.0
mkdir -p "$DB_DIR" && cd "$DB_DIR"

echo "[1/4] geNomad database..."
[ -d "$GENOMAD_DB" ] || genomad download-database "$DB_DIR"

echo "[2/4] CheckV database v1.5 (+ rebuild with pinned DIAMOND 2.0.15)..."
if [ ! -d "$CHECKV_DB" ]; then
  checkv download_database "$DB_DIR"
  ( cd "${CHECKV_DB}/genome_db" && diamond makedb --in checkv_reps.faa --db checkv_reps --threads 4 )
fi

echo "[3/4] DefenseFinder models 3.1.0 + CasFinder 3.1.0..."
if [ ! -d "$DF_MODELS" ]; then
  mkdir -p "$DF_MODELS" && cd "$DF_MODELS"
  curl -L -o models.tar.gz \
    https://github.com/mdmparis/defense-finder-models/releases/download/3.1.0/defense-finder-models-3.1.0.tar.gz
  tar xzf models.tar.gz && rm -f models.tar.gz
  git clone --branch 3.1.0 --depth 1 https://github.com/macsy-models/CasFinder.git
  grep -q '^vers:' CasFinder/metadata.yml || echo "vers: '3.1.0'" >> CasFinder/metadata.yml
  cd "$DB_DIR"
fi

echo "[4/4] GTDB r232 bacterial reference tree (bac120)..."
mkdir -p "$(dirname "$GTDB_TREE")"
if [ ! -s "$GTDB_TREE" ]; then
  curl -L -o "${GTDB_TREE}.gz" https://data.gtdb.ecogenomic.org/releases/release232/232.0/bac120_r232.tree.gz
  gunzip -f "${GTDB_TREE}.gz"
fi
echo "DONE. Databases in $DB_DIR"; du -sh "$DB_DIR"/* 2>/dev/null
