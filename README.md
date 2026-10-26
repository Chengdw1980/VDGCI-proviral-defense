# VDGCI pipeline and analysis code (v3.2)

Code accompanying **"Defense enrichment, phylogenetic structure and spatial
organization of proviral DNA across bacterial orders and classes"**

Davis W. Cheng, Ph.D. · Rocky Institute of Biotechnology, Clovis, California, USA
ORCID [0009-0001-1306-8923](https://orcid.org/0009-0001-1306-8923)

Companion data deposit: aggregated validation-cohort outputs (see *Data
availability* in the article). Genome sequences are public at NCBI under the
accessions listed in `config/`.

---

## Scope

The pipeline quantifies integrated viral sequence and antiphage defense content
across bacterial genomes, and localizes each defense system relative to
predicted proviral boundaries. It was applied to two cohorts:

| Cohort | Genomes | Purpose |
|---|---|---|
| Discovery | 48 complete RefSeq assemblies | hypothesis generation |
| Validation | 5,054 GTDB r232 family representatives | prespecified tests V1–V8 and exploratory taxon-resolved analyses |

## Software versions

Results in the article were produced with geNomad 1.12.0 (database v1.9),
CheckV, PhiSpy, DefenseFinder 3.0.0 (models 3.1.0; CasFinder 3.1.0), CheckM2,
pyrodigal, and GTDB release 232. Conda specifications are in `environment.yml`,
`environment_df.yml` and `environment_diamond.yml`.

Computation was performed on Purdue Anvil under ACCESS allocation BIO260593.
`docs/README_ANVIL.md` documents the cluster-specific execution path and
`docs/INSTALL_NOTES.md` the database setup.

## Layout

```
config/        accession lists, cohort manifests with full GTDB taxonomy,
               mechanism-class lookup, pipeline configuration
scripts/       01 database setup .. 07 cost report — the per-genome pipeline
               and the aggregation/validation statistics (05, 06)
analysis/      post-pipeline analyses: defense-system mapping and positions,
               robustness and audit statistics, taxon-resolved enrichment,
               phylogenetic signal, power
figures/       figure-generating scripts
anvil_jobs/    SLURM submission scripts used on Anvil
reference/     curated RefSeq bacteriophage accession set
docs/          cluster execution notes, installation notes, pipeline README v2
```

## Reproducing the published analyses

1. `scripts/01_setup_databases.sh`, `scripts/02_fetch_genomes.sh` — databases and genomes.
2. `scripts/03_*` — per-genome geNomad, DefenseFinder, PhiSpy and CheckV.
3. `scripts/04_gtdb_tree.py` — prune the GTDB bac120 tree to the cohort.
4. `scripts/05_analyze.py`, `scripts/06_extended_analysis.py` — cohort aggregation and the prespecified validation tests.
5. `analysis/defmap.py` → `defstats.py` → `defmodel.py` → `defaudit.py` → `defrobust.py` — defense-system coordinates, positional statistics and the audit of system membership.
6. `analysis/defclass.py` — class- and order-resolved enrichment, class × system-type associations, class-resolved positional tests and provirus–tRNA adjacency (Table 2, Fig. 7, Supplementary Table S13a–S13d).
7. `figures/` — display items.

## Version note (important)

**v3.2** corrects a defect in which DefenseFinder emitted exact duplicate
system records, inflating system counts. `analysis/patch_dedup.py` documents
the correction applied to `scripts/05_analyze.py` and
`scripts/06_extended_analysis.py`; those scripts in this release already
contain it. Every defense analysis reported in the article was recomputed after
the correction. Results produced with earlier versions are superseded.

## Licence

MIT (see `LICENSE`). Please cite the article and this software; see
`CITATION.cff`.
