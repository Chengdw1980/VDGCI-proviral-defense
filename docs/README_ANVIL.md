# Running the VDGCI pipeline on Purdue Anvil (v3, 23 October 2026)

**New in v3:** (1) each genome also keeps geNomad's plasmid calls and per-gene viral annotations (small files; needed for the plasmid comparison and the prophage-decay analysis); (2) new `scripts/06_extended_analysis.py` runs the pre-specified statistics for the paper (mechanism-class enrichment with genome-level bootstrap CIs, per-phylum replication, defense density by CheckV completeness tier, phylogenetic signal on the GTDB tree, genome-level model, provirus-vs-plasmid comparison); (3) `config/defense_mechanism_map.tsv` assigns each defense type to a mechanism class and can be edited.

Project: ACCESS **BIO260593** — PI Davis W. Cheng, Rocky Institute of Biotechnology  
Resource: **Purdue Anvil CPU, 50,000 SUs** (1 SU = 1 core-hour)

Copy each grey command box exactly. Lines starting with `#` are comments.


## 0. Before you start
* Wait for the email saying your exchange is approved and your **Anvil account** exists.
* Your Anvil username is **x-dcheng3**.

## 1. Open a terminal on Anvil (no software needed on your PC)
1. Go to **https://ondemand.anvil.rcac.purdue.edu** and log in with your ACCESS ID (dcheng3) + Duo.
2. Top menu **Clusters → Anvil Shell Access**. A black terminal window opens.

## 2. Upload this pipeline
1. In OnDemand, top menu **Files → Home Directory → Upload**, choose `VDGCI_pipeline_Anvil_v2.zip`.
2. In the terminal:
```bash
cd $PROJECT
mv ~/VDGCI_pipeline_Anvil_v2.zip . && unzip -q VDGCI_pipeline_Anvil_v2.zip
cd VDGCI_cluster_pipeline
```

## 3. Tell the pipeline your allocation name
```bash
mybalance
```
Note the CPU account name it prints (it ends in `-cpu`). Then (replace NAME):
```bash
sed -i 's/CHANGE_ME-cpu/NAME/' config/config.sh
```

## 4. Install the software (once, ~30–60 min)
```bash
# Anvil's conda adds the 'defaults' channel, which makes the solve fail. Use only conda-forge + bioconda:
printf 'channels:\n  - conda-forge\n  - bioconda\nchannel_priority: strict\n' > ~/.condarc
nohup bash scripts/00_install_envs.sh > $PROJECT/install_env.log 2>&1 &   # builds 3 environments + wrappers (v3.1)
tail -f $PROJECT/install_env.log     # wait for: ALL ENVIRONMENTS READY  (Ctrl+C stops watching, not the install)
```

## 5. Download the databases (once, ~1 h, ~15 GB)
```bash
module load anaconda; eval "$(conda shell.bash hook)"; conda activate $PROJECT/vdgci/env
bash scripts/01_setup_databases.sh
```

## 6. TEST RUN on 10 control genomes (always do this first)
```bash
export ACCESSION_LIST=$PWD/config/accessions_test10.txt
bash scripts/02_fetch_genomes.sh
mkdir -p logs_slurm
sbatch -A NAME scripts/03_test_shared.sbatch
squeue -u $USER          # repeat until the jobs disappear (usually 1–3 h)
source config/config.sh
python3 scripts/05_analyze.py
python3 scripts/07_cost_report.py 5054
```
**Check the controls before going further:** *E. coli* K-12 (GCF_000005845.2) must show
proviruses; *Mycoplasma genitalium* (GCF_000027325.1) must show zero. The cost report
projects how many SUs the full run will use — send it to Claude to review.

## 7. FULL RUN — 5,054 genomes (only after the test passes)
```bash
export ACCESSION_LIST=$PWD/config/accessions_HQ_5054.txt
bash scripts/02_fetch_genomes.sh            # ~1–3 h on the login node
bash scripts/03_submit_chunks.sh            # submits ~16 whole-node jobs
squeue -u $USER                             # watch progress
ls $PROJECT/vdgci/HQ_5054/done | wc -l   # genomes finished
```
If anything stops early, just run `bash scripts/03_submit_chunks.sh` again — it only
resubmits unfinished genomes.

## 8. Analysis and host tree
```bash
source config/config.sh
python3 scripts/04_gtdb_tree.py             # GTDB r232 tree pruned to the cohort
python3 scripts/05_analyze.py
python3 scripts/06_extended_analysis.py     # pre-specified statistics (about 30-90 min on a login/compute node)
```
`06` writes to `$PROJECT/vdgci/HQ_5054/results/extended/`. Run it inside a short interactive job if the login node is busy:
```bash
sinteractive -A NAME -p shared -n 8 --mem=32G -t 2:00:00
```
Results land in `$PROJECT/vdgci/HQ_5054/results/` (safe storage). Download them with
OnDemand **Files**.


## What changed from v1 (and why)
| Issue at genome scale | v1 behaviour | v2 fix |
|---|---|---|
| 83% of the cohort are GenBank drafts/MAGs; most lack a 16S gene | 16S tree would drop them | Host tree = **GTDB r232 bac120 tree** pruned to the cohort (5,054/5,054 placed; tested) |
| PhiSpy needs real NCBI annotation (only 2,562/5,054 = 50.7% have it) | consensus % diluted by genomes PhiSpy never ran on | consensus computed **only on annotated genomes**; reported with its n |
| Fragmented MAGs split prophages (829 genomes have contig N50 < 20 kb) | not controlled | contigs + N50 recorded per genome; enrichment repeated with **N50 ≥ 20 kb** |
| Anvil job limits (~128 queued jobs per partition) | 5,000-task array | **16 genomes per 128-core node**, ~16 node-jobs total |
| Anvil scratch is purged after 30 days | results on scratch | small result files copied to **$PROJECT** as each genome finishes |
| Anvil charges by memory share (~2 GB/core) | 8 cores + 32 GB = billed as 16 | 8 cores + 16 GB, geNomad `--splits 8` = billed as 8 |
| 5,054 single NCBI downloads | slow, fragile | bulk dehydrated download in batches of 500 |

Tested in Claude's workspace: GTDB pruning (5,054/5,054), analysis script on synthetic
data with hand-checked answers, and PhiSpy on *E. coli* K-12 (4 prophages with NCBI
annotation vs 0 from bare gene calls — the reason PhiSpy is limited to annotated genomes).
Not yet tested: the tools themselves on Anvil — that is what step 6 is for.
