# Install notes and known gotchas

## 0. (v3.1, 23 Oct 2026) Three environments, not one
The single environment.yml could not be solved on Anvil: DIAMOND 2.0.15 conflicts
with geNomad's TensorFlow dependencies, and Anvil's default `defaults` channel adds
more conflicts. Fix: (a) put `channels: [conda-forge, bioconda]` and
`channel_priority: strict` in ~/.condarc; (b) run `scripts/00_install_envs.sh`,
which builds `env` (main), `env_diamond` (DIAMOND 2.0.15 only) and `env_df`
(DefenseFinder 3.0.0 = PyPI package `mdmparis-defense-finder`, MacSyFinder 2.1.4,
HMMER) and writes wrappers to `$TOOLS_BIN` that scripts 01 and 03 put first on PATH.
The PyPI name `defense-finder` used in earlier notes does not exist.
DefenseFinder 3.0.0 + models 3.1.0 + CasFinder 3.1.0 was re-tested on E. coli K-12
proteins (finds type I-E Cas, RM type I, McrBC, MazEF, Hachiman, Lit, RnlAB).

These are the exact traps we hit building the study pipeline. Reproducing them
on your cluster will save hours.

## 1. DefenseFinder models (the one most likely to bite you)
DefenseFinder needs two model sets: `defense-finder-models` (v3.1.0) and
`CasFinder` (v3.1.0). On a locked-down/proxied network the automatic
`defense-finder update` fails because it pulls from `api.github.com/codeload`,
which is frequently blocked. Work around it:

```bash
mkdir -p db/df_models && cd db/df_models
# defense-finder models 3.1.0 (release asset, not codeload):
curl -L -o models.tar.gz \
  https://github.com/mdmparis/defense-finder-models/releases/download/3.1.0/defense-finder-models-3.1.0.tar.gz
tar xzf models.tar.gz
# CasFinder 3.1.0 — you MUST use the tag, not the default branch:
git clone --branch 3.1.0 --depth 1 https://github.com/macsy-models/CasFinder.git
```
The CasFinder **default branch uses MacSyFinder grammar 2.1, which is
incompatible** with the DefenseFinder version used here; **tag 3.1.0 uses
grammar 2.0** and works. If MacSyFinder complains about a missing `vers:` field,
add `vers: '3.1.0'` to `CasFinder/metadata.yml`.

Then always run DefenseFinder with:
`defense-finder run --skip-model-version-check --models-dir db/df_models ...`
(`--skip-model-version-check` avoids an online SSL version check that fails
offline).

## 2. CheckV + DIAMOND version
CheckV calls `diamond ... -k 10000`. **DIAMOND 2.1.9 segfaults on that flag.**
Pin **DIAMOND 2.0.15** (done in `environment.yml`). After downloading the CheckV
DB you must (re)build the reps database with the pinned DIAMOND:

```bash
cd db/checkv-db-v1.5/genome_db
diamond makedb --in checkv_reps.faa --db checkv_reps --threads $THREADS
```

## 3. geNomad memory (OOM on marker search)
`genomad end-to-end` runs an mmseqs2 search against ~228k profiles. On low-RAM
nodes this dies with "Can not allocate entries memory". Fix by splitting:
`--splits N` (we used 4–8). More splits = less RAM, slightly slower. Set via
`GENOMAD_SPLITS` in `config/config.sh`. On a real cluster node (>=32 GB) you can
usually use `--splits 2` or omit it.

## 4. Sequence naming for batched runs
If you batch many genomes into one geNomad run to amortize the marker search,
rename every contig `>{ACCESSION}__{original_id}` so results can be mapped back
to genomes (double underscore is the delimiter the analysis script splits on).
The SLURM array approach in this package runs **one genome per task instead**,
which is simpler and scales linearly — batching is only worth it on a single
shared-memory machine.

## 5. Data provenance
Genomes: NCBI RefSeq assemblies (GCF_*), retrieved with `datasets`.
Clean phage reference: RefSeq viral release, classified to host domain via the
Virus-Host DB (`reference/refseq_phage_accessions.csv`, 6,649 phage/archaeal-virus
genomes). Do **not** use an unfiltered "all viruses" FASTA — ours was 96% one
eukaryotic virus (BVDV) before filtering.
