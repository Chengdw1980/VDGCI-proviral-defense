#!/bin/bash
# Run every per-genome tool for ONE accession, then copy the small result files
# to $KEEP (project storage) so Anvil's 30-day scratch purge can never lose them.
# Idempotent: finished steps are skipped, so a re-run resumes where it stopped.
#   Usage: scripts/03_run_per_genome.sh <ACCESSION>
set -uo pipefail
source "$(dirname "$0")/../config/config.sh"
export PATH="${TOOLS_BIN}:${PATH}"   # pinned DIAMOND 2.0.15 + DefenseFinder 3.0.0
ACC="$1"; T0=$(date +%s)
[ -s "${KEEP}/done/${ACC}" ] && { echo "$ACC already done"; exit 0; }
FNA="${GENOMES_DIR}/${ACC}.fna"; GBFF="${GBFF_DIR}/${ACC}.gbff"
[ -s "$FNA" ] || { echo "$ACC: no genome FASTA (fetch failed?)"; echo "$ACC" >> "${KEEP}/failed_no_fasta.txt"; exit 0; }
echo "######## $ACC start $(date -Is) on $(hostname)"

# 1. geNomad — integrated viral elements (+ its own gene calls, used for all genomes)
GOUT="${GENOMAD_OUT}/${ACC}"
if [ ! -s "${GOUT}/${ACC}_summary/${ACC}_virus_summary.tsv" ]; then
  genomad end-to-end --cleanup --splits "$GENOMAD_SPLITS" --threads "$THREADS" \
      "$FNA" "$GOUT" "$GENOMAD_DB" > "${LOGS}/${ACC}.genomad.log" 2>&1 \
      || echo "$ACC geNomad error (see log)"
fi
# 2. DefenseFinder — on geNomad's proteins, so gene coordinates match provirus coordinates
GPROT="${GOUT}/${ACC}_annotate/${ACC}_proteins.faa"
DOUT="${DF_OUT}/${ACC}"
if [ -s "$GPROT" ] && [ ! -s "${DOUT}/${ACC}_proteins_defense_finder_systems.tsv" ] \
                   && [ ! -s "${DOUT}/${ACC}_defense_finder_systems.tsv" ]; then
  defense-finder run --skip-model-version-check --models-dir "$DF_MODELS" \
      -w "$THREADS" -o "$DOUT" "$GPROT" > "${LOGS}/${ACC}.df.log" 2>&1 || echo "$ACC DefenseFinder error"
fi
# 3. PhiSpy — independent caller; needs a real NCBI annotation, so it runs only
#    where a GBFF exists (it relies partly on gene-function names; tested: it
#    finds 0 prophages in E. coli K-12 from bare gene calls vs 4 with NCBI's).
POUT="${PHISPY_OUT}/${ACC}"
if [ -s "$GBFF" ] && grep -q "^     CDS " "$GBFF" && [ ! -e "${POUT}/prophage_coordinates.tsv" ]; then   # v3.1: only annotated GBFFs (about 1/3 of downloaded GBFFs have no CDS)
  PhiSpy.py "$GBFF" -o "$POUT" --threads "$THREADS" > "${LOGS}/${ACC}.phispy.log" 2>&1 \
     || { mkdir -p "$POUT"; touch "${POUT}/prophage_coordinates.tsv"; }
  [ -e "${POUT}/prophage_coordinates.tsv" ] || touch "${POUT}/prophage_coordinates.tsv"
fi
# 4. CheckV — completeness of this genome's proviruses
PROV="${GOUT}/${ACC}_summary/${ACC}_virus.fna"; COUT="${CHECKV_OUT}/${ACC}"
if [ -s "$PROV" ] && [ ! -s "${COUT}/quality_summary.tsv" ]; then
  checkv end_to_end "$PROV" "$COUT" -d "$CHECKV_DB" -t "$THREADS" > "${LOGS}/${ACC}.checkv.log" 2>&1 \
     || echo "$ACC CheckV error"
fi

# 5. copy the small, analysis-ready files to safe project storage
K="${KEEP}/per_genome/${ACC}"; mkdir -p "$K" "${KEEP}/done"
cp -f "${GOUT}/${ACC}_summary/${ACC}_virus_summary.tsv" "$K/" 2>/dev/null
cp -f "$PROV" "$K/" 2>/dev/null
# v3: keep geNomad plasmid calls (MGE comparison class) and per-gene viral annotations
#     (structural/hallmark genes for the decay and satellite analyses)
cp -f "${GOUT}/${ACC}_summary/${ACC}_plasmid_summary.tsv" "${GOUT}/${ACC}_summary/${ACC}_virus_genes.tsv" "$K/" 2>/dev/null
grep '^>' "$GPROT" > "$K/${ACC}_protein_headers.txt" 2>/dev/null      # gene coordinates
cp -f "${DOUT}"/*defense_finder_systems.tsv "${DOUT}"/*defense_finder_genes.tsv "$K/" 2>/dev/null
[ -e "${POUT}/prophage_coordinates.tsv" ] && cp -f "${POUT}/prophage_coordinates.tsv" "$K/phispy_prophage_coordinates.tsv"
cp -f "${COUT}/quality_summary.tsv" "$K/checkv_quality_summary.tsv" 2>/dev/null
# genome size and fragmentation (contigs, N50) — covariates for the analysis
python3 - "$FNA" > "$K/genome_stats.tsv" <<'PY'
import sys
L=[];n=0
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if n: L.append(n)
        n=0
    else: n+=len(line.strip())
if n: L.append(n)
L.sort(reverse=True); tot=sum(L); c=0; n50=0
for x in L:
    c+=x
    if c>=tot/2: n50=x; break
print("genome_bp\tn_contigs\tcontig_n50"); print(f"{tot}\t{len(L)}\t{n50}")
PY
T1=$(date +%s)
# mark done only if geNomad produced its summary (the core result)
if [ -s "$K/${ACC}_virus_summary.tsv" ]; then
  echo -e "${ACC}\t$((T1-T0))\t${THREADS}" > "${KEEP}/done/${ACC}"
  echo "######## $ACC done in $((T1-T0)) s"
else
  echo "$ACC" >> "${KEEP}/failed_genomad.txt"; echo "######## $ACC FAILED (no geNomad summary)"
fi
