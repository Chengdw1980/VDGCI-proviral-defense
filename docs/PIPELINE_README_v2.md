> **Version 2 (October 2026) — Anvil-ready.** For step-by-step instructions on
> Purdue Anvil (ACCESS project BIO260593) read **README_ANVIL.md** first. It
> supersedes the SLURM-array instructions below, which describe v1.

# Viral-Derived Genomic Contribution (VDGCI) pipeline

Reproducible, cluster-ready pipeline for the study *"Viral-derived DNA is a
labile, defense-enriched compartment of bacterial genomes."* It quantifies the
integrated viral compartment of bacterial genomes and its anti-phage defense
content, and scales from the 48-genome published cohort to **thousands of
genomes** via a SLURM job array (one genome per task).

Author: Davis W. Cheng, Ph.D. — Rocky Institute of Biotechnology.

---

## What it does (per genome, in parallel)
1. **geNomad** — detect integrated viral elements (proviruses) and plasmids.
2. **DefenseFinder** — annotate anti-phage defense systems (run on geNomad's
   protein calls so each system can be placed inside/outside a provirus).
3. **PhiSpy** — independent, composition-based prophage caller (corroboration).
4. **CheckV** — completeness/quality of each provirus.

Then, once (single jobs):
5. **16S phylogeny** — MAFFT + FastTree, for phylogenetic-signal tests.
6. **Aggregation** — per-genome VDGCI, cohort defense-enrichment binomial test
   and gene density, geNomad↔PhiSpy consensus, CheckV spectrum, output tables.

## Layout
```
config/     config.sh (edit paths/params), accessions.txt, cohort_manifest.tsv
scripts/    01_setup_databases.sh  02_fetch_genomes.sh
            03_run_per_genome.sh   03_array.sbatch   (the parallel core)
            04_phylogeny.sh        05_analyze.py
reference/  refseq_phage_accessions.csv (clean 6,649-genome phage reference set)
environment.yml   INSTALL_NOTES.md   (READ THIS — the version traps are real)
```

## Quick start (cluster)
```bash
# 0. environment
mamba env create -f environment.yml && conda activate vdgci

# 1. choose your genomes: replace config/accessions.txt with your GCF_* list.
#    Start point provided: config/accessions_48.txt (the published cohort).
cp config/accessions_48.txt config/accessions.txt      # or your big list

# 2. one-time database download (login/transfer node, needs network)
bash scripts/01_setup_databases.sh

# 3. fetch genomes (login/transfer node, needs network)
bash scripts/02_fetch_genomes.sh

# 4. run the per-genome array — the scaling step
source config/config.sh
N=$(grep -cve '^\s*$' config/accessions.txt)
sbatch --array=1-${N}%100 scripts/03_array.sbatch      # %100 = up to 100 concurrent genomes

# 5. after the array finishes: phylogeny + aggregation (single jobs)
bash scripts/04_phylogeny.sh
python3 scripts/05_analyze.py     # writes tables to $WORK/results/
```

## Scaling notes
- **The array is the whole story.** Each genome is independent, so wall-clock is
  ~(one genome's time) regardless of cohort size, bounded by your queue/`%`
  concurrency. 5,000 genomes at 100 concurrent ≈ 50 waves.
- Tune `--mem`, `--cpus-per-task`, `--time` in `03_array.sbatch` and `THREADS`,
  `GENOMAD_SPLITS` in `config/config.sh`. geNomad's marker search is the memory
  driver; on ≥32 GB nodes `GENOMAD_SPLITS=2` is usually fine.
- Everything is **idempotent** — resubmit the array to resume after failures;
  finished genomes are skipped.
- Put `$WORK` on fast scratch, not home. Expect ~50–200 MB of output per genome.

## Reproducing the 48-genome result
Use `config/accessions_48.txt` as `config/accessions.txt`. The aggregator prints
the headline numbers (viral-derived %, defense enrichment fold + binomial p,
consensus %, CheckV tiers) and writes the same tables used in the manuscript.
Small numeric differences from the published values are expected if tool/DB
versions differ from those pinned in `environment.yml`.

## Validation & honesty note
This pipeline computes correlational, annotation-based quantities. Defense-system
"functions" are **predicted/annotated**, not experimentally demonstrated. The
VDGCI is a descriptive composite index, not a validated evolutionary rate. Treat
outputs accordingly in any manuscript claims.
