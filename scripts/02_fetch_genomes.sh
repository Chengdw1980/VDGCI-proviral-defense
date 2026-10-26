#!/bin/bash
# Download every genome in $ACCESSION_LIST from NCBI in bulk (dehydrated
# package + rehydrate = far faster and gentler on NCBI than 5,000 single calls).
# Gets nucleotide FASTA for all, plus protein FASTA and GBFF where NCBI has an
# annotation (~51% of the HQ cohort; GenBank-only MAGs often have none).
# Run on an Anvil LOGIN node (needs internet). Safe to re-run: skips finished.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
DL="${WORK}/_ncbi_bulk"; mkdir -p "$DL"
TODO="${DL}/todo.txt"; : > "$TODO"
while read -r ACC; do
  [ -z "$ACC" ] && continue
  [ -s "${GENOMES_DIR}/${ACC}.fna" ] || echo "$ACC" >> "$TODO"
done < "$ACCESSION_LIST"
echo "Genomes still to fetch: $(wc -l < "$TODO") of ${N_GENOMES}"
[ -s "$TODO" ] || { echo "Nothing to do."; exit 0; }

# work in batches of 500 so one failure never loses the whole download
split -l 500 -d "$TODO" "${DL}/batch_"
for B in "${DL}"/batch_??; do
  echo "=== $(basename "$B"): $(wc -l < "$B") genomes ==="
  rm -rf "${B}_pkg" "${B}.zip"
  datasets download genome accession --inputfile "$B" \
      --include genome,protein,gbff --dehydrated --filename "${B}.zip"
  unzip -q -o "${B}.zip" -d "${B}_pkg"
  datasets rehydrate --directory "${B}_pkg" --max-workers 8
  for D in "${B}_pkg"/ncbi_dataset/data/GC*; do
    ACC=$(basename "$D")
    cat "$D"/*_genomic.fna > "${GENOMES_DIR}/${ACC}.fna" 2>/dev/null || true
    [ -s "$D/protein.faa" ]    && cp "$D/protein.faa"    "${PROT_DIR}/${ACC}.faa"
    [ -s "$D/genomic.gbff" ]   && cp "$D/genomic.gbff"   "${GBFF_DIR}/${ACC}.gbff"
  done
  rm -rf "${B}_pkg" "${B}.zip" "$B"
done

n_fna=$(ls "${GENOMES_DIR}" | wc -l); n_gb=$(ls "${GBFF_DIR}" | wc -l)
echo "Done. nucleotide: ${n_fna}/${N_GENOMES}   NCBI-annotated (GBFF): ${n_gb}"
# list anything still missing so it can be retried
: > "${KEEP}/fetch_missing.txt"
while read -r ACC; do [ -z "$ACC" ] || [ -s "${GENOMES_DIR}/${ACC}.fna" ] || echo "$ACC" >> "${KEEP}/fetch_missing.txt"; done < "$ACCESSION_LIST"
echo "Missing (see ${KEEP}/fetch_missing.txt): $(wc -l < "${KEEP}/fetch_missing.txt")"
