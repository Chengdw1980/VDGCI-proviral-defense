#!/usr/bin/env python3
"""
Aggregate per-genome results (from $KEEP/per_genome, i.e. safe project storage)
into the study tables and cohort statistics:
  - per-genome viral-derived fraction, element count, hallmarks, VDGCI
  - genome size and assembly fragmentation (contigs, N50) as covariates
  - defense-in-viral-DNA enrichment (binomial test + gene density), computed on
    genomes with BOTH geNomad and DefenseFinder results; repeated excluding
    fragmented assemblies (contig N50 < 20 kb) as a sensitivity analysis
  - geNomad<->PhiSpy corroboration, computed ONLY on genomes where PhiSpy ran
    (NCBI-annotated genomes), so the rate is not diluted by unannotated MAGs
  - CheckV completeness tiers
Usage:  source config/config.sh && python3 scripts/05_analyze.py
"""
import os, csv, glob, json, sys
from collections import defaultdict, Counter
from scipy import stats

def env(k, d=None):
    v = os.environ.get(k, d)
    if v is None: sys.exit(f"Missing {k}; run `source config/config.sh` first.")
    return v
KEEP = env("KEEP"); RESULTS = env("RESULTS"); ROOT = env("PIPELINE_ROOT")
OVERLAP = float(os.environ.get("CONSENSUS_OVERLAP", "0.5"))
N50_MIN = int(os.environ.get("SENS_MIN_N50", "20000"))
PG = os.path.join(KEEP, "per_genome"); os.makedirs(RESULTS, exist_ok=True)

# ---------- taxonomy / metadata ----------
meta = {}
for cand in [os.environ.get("MANIFEST", ""), os.path.join(ROOT, "config", "cohort_manifest.tsv"),
             os.path.join(ROOT, "config", "cohort_manifest_full.tsv"),
             os.path.join(ROOT, "config", "cohort_manifest_48.tsv")]:
    if cand and os.path.exists(cand):
        for r in csv.DictReader(open(cand), delimiter="\t"):
            meta.setdefault(r["accession"], {"species": r.get("species", ""), "phylum": r.get("phylum", ""),
                "class": r.get("class", ""), "family": r.get("family", ""), "lifestyle": r.get("lifestyle", "")})
blank = {"species": "", "phylum": "", "class": "", "family": "", "lifestyle": ""}

genomes = {}; prov = defaultdict(list); gene_xy = {}; prov_rows = []
for d in sorted(glob.glob(os.path.join(PG, "*"))):
    acc = os.path.basename(d)
    vs = os.path.join(d, f"{acc}_virus_summary.tsv")
    if not os.path.exists(vs): continue
    gs = {"genome_bp": 0, "n_contigs": 0, "contig_n50": 0}
    gsf = os.path.join(d, "genome_stats.tsv")
    if os.path.exists(gsf):
        r = next(csv.DictReader(open(gsf), delimiter="\t")); gs = {k: int(r[k]) for k in gs}
    g = dict(acc=acc, viral_bp=0, n_prov=0, hall=0, **gs, **meta.get(acc, blank))
    for r in csv.DictReader(open(vs), delimiter="\t"):
        sn = r["seq_name"]; L = int(r["length"]); co = (r.get("coordinates") or "").strip()
        hall = int(float(r.get("n_hallmarks", 0) or 0))
        if co and co.upper() != "NA" and "-" in co:
            contig = sn.split("|")[0]; s, e = map(int, co.split("-"))
        else:
            contig, s, e = sn, 1, L
        prov[(acc, contig)].append((s, e)); g["viral_bp"] += e - s + 1; g["n_prov"] += 1; g["hall"] += hall
        prov_rows.append([acc, g["species"], g["phylum"], contig, s, e, e - s + 1, hall,
                          r.get("virus_score", ""), r.get("taxonomy", "")])
    hf = os.path.join(d, f"{acc}_protein_headers.txt")
    if os.path.exists(hf):
        for line in open(hf):
            p = line[1:].split("#")
            try: gene_xy[(acc, p[0].strip())] = (p[0].strip().rsplit("_", 1)[0], int(p[1]), int(p[2]))
            except (IndexError, ValueError): pass
    g["pct_viral"] = 100 * g["viral_bp"] / g["genome_bp"] if g["genome_bp"] else 0
    genomes[acc] = g
if not genomes: sys.exit(f"No per-genome results in {PG}")
accs = sorted(genomes)

def norm(v):
    lo, hi = min(v), max(v); return [0.0 if hi == lo else (x - lo) / (hi - lo) for x in v]
for a, x, y, z in zip(accs, norm([genomes[a]["pct_viral"] for a in accs]),
                      norm([genomes[a]["hall"] for a in accs]), norm([genomes[a]["n_prov"] for a in accs])):
    genomes[a]["VDGCI"] = round((x + y + z) / 3, 4)

def in_prov(acc, contig, s, e):
    return any(s <= pe and e >= ps for ps, pe in prov.get((acc, contig), []))

# ---------- DefenseFinder ----------
sys_rec = []; per_def = {}
for a in accs:
    d = os.path.join(PG, a)
    sf = glob.glob(os.path.join(d, "*defense_finder_systems.tsv")); gf = glob.glob(os.path.join(d, "*defense_finder_genes.tsv"))
    if not sf or not gf: continue
    sg = defaultdict(list)
    for r in csv.DictReader(open(gf[0]), delimiter="\t"): sg[r.get("sys_id", "")].append(r.get("hit_id", ""))
    sg = {k: list(dict.fromkeys(v)) for k, v in sg.items()}; seen = set()  # de-duplicate (v3.2)
    c = {"n_sys": 0, "n_sys_prov": 0, "n_gene": 0, "n_gene_prov": 0}
    for r in csv.DictReader(open(sf[0]), delimiter="\t"):
        if r.get("sys_id", "") in seen: continue  # DefenseFinder can list a system twice
        seen.add(r.get("sys_id", ""))
        inp = False
        for pid in sg.get(r.get("sys_id", ""), []):
            xy = gene_xy.get((a, pid))
            if xy:
                c["n_gene"] += 1
                if in_prov(a, *xy): c["n_gene_prov"] += 1; inp = True
        c["n_sys"] += 1; c["n_sys_prov"] += inp
        sys_rec.append((a, r.get("sys_id", ""), r.get("type", ""), r.get("subtype", ""), inp))
    per_def[a] = c

def enrichment(subset):
    sub = [a for a in subset if a in per_def]
    S = sum(per_def[a]["n_sys"] for a in sub); K = sum(per_def[a]["n_sys_prov"] for a in sub)
    bp = sum(genomes[a]["genome_bp"] for a in sub); vbp = sum(genomes[a]["viral_bp"] for a in sub)
    if not S or not vbp or not bp: return {"n_genomes": len(sub), "S": S}
    p0 = vbp / bp; bt = stats.binomtest(K, S, p0, alternative="greater")
    G = sum(per_def[a]["n_gene"] for a in sub); Gp = sum(per_def[a]["n_gene_prov"] for a in sub)
    dp = Gp / (vbp / 1e6); dr = (G - Gp) / ((bp - vbp) / 1e6)
    return dict(n_genomes=len(sub), S=S, K=K, pct_systems_in_viral=100 * K / S, pct_sequence_viral=100 * p0,
                fold=(K / S) / p0, binom_p=bt.pvalue, density_viral_per_Mb=dp, density_rest_per_Mb=dr,
                density_ratio=dp / dr if dr else None)
enr_all = enrichment(accs)
enr_sens = enrichment([a for a in accs if genomes[a]["contig_n50"] >= N50_MIN])

# ---------- PhiSpy corroboration (only where PhiSpy ran) ----------
def ovf(a, b):
    ov = max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1); m = min(a[1] - a[0] + 1, b[1] - b[0] + 1)
    return ov / m if m > 0 else 0
ran = {}; 
for a in accs:
    pc = os.path.join(PG, a, "phispy_prophage_coordinates.tsv")
    if os.path.exists(pc):
        ran[a] = []
        for line in open(pc):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                try: ran[a].append((p[1], int(p[2]), int(p[3])))
                except ValueError: pass
gt = gc = 0
for (a, contig), ivs in prov.items():
    if a not in ran: continue
    for iv in ivs:
        gt += 1; gc += any(c == contig and ovf(iv, (s, e)) >= OVERLAP for c, s, e in ran[a])
consensus = dict(genomes_with_phispy=len(ran), genomad_proviruses_in_those=gt, corroborated=gc,
                 pct=100 * gc / gt if gt else None, phispy_total=sum(len(v) for v in ran.values()))

# ---------- CheckV ----------
cv = Counter(); cv_rows = []
for a in accs:
    q = os.path.join(PG, a, "checkv_quality_summary.tsv")
    if os.path.exists(q):
        for r in csv.DictReader(open(q), delimiter="\t"):
            cv[r.get("checkv_quality", "NA")] += 1
            cv_rows.append([a, r.get("contig_id", ""), r.get("completeness", ""), r.get("checkv_quality", "")])

# ---------- write ----------
def w(name, header, rows, delim=","):
    with open(os.path.join(RESULTS, name), "w", newline="") as f:
        x = csv.writer(f, delimiter=delim); x.writerow(header); x.writerows(rows)
w("cohort_VDGCI.csv", ["accession", "species", "phylum", "class", "family", "lifestyle", "genome_bp", "n_contigs",
  "contig_n50", "viral_bp", "pct_viral", "n_viral_elements", "n_hallmarks", "VDGCI", "phispy_ran"],
  [[a, g["species"], g["phylum"], g["class"], g["family"], g["lifestyle"], g["genome_bp"], g["n_contigs"], g["contig_n50"],
    g["viral_bp"], round(g["pct_viral"], 4), g["n_prov"], g["hall"], g["VDGCI"], a in ran]
   for a, g in sorted(genomes.items(), key=lambda kv: -kv[1]["VDGCI"])])
w("provirus_detail.csv", ["accession", "species", "phylum", "contig", "start", "end", "length_bp", "n_hallmarks",
  "virus_score", "taxonomy"], prov_rows)
w("defense_per_genome.csv", ["accession", "n_defense_systems", "n_in_viral_dna", "n_defense_genes", "n_genes_in_viral_dna"],
  [[a, c["n_sys"], c["n_sys_prov"], c["n_gene"], c["n_gene_prov"]] for a, c in sorted(per_def.items())])
w("defense_systems_detail.csv", ["accession", "sys_id", "type", "subtype", "in_viral_dna"], sys_rec)
w("checkv_summary.tsv", ["accession", "contig_id", "completeness", "checkv_quality"], cv_rows, "\t")
by_type = defaultdict(lambda: [0, 0])
for _, _, t, _, inp in sys_rec: by_type[t][0] += 1; by_type[t][1] += inp
w("defense_by_type.csv", ["type", "n_systems", "n_in_viral_dna", "pct_in_viral_dna"],
  sorted([[t, n, k, round(100 * k / n, 2)] for t, (n, k) in by_type.items()], key=lambda r: -r[1]))
total_bp = sum(g["genome_bp"] for g in genomes.values()); vbp = sum(g["viral_bp"] for g in genomes.values())
summary = dict(n_genomes=len(accs), n_viral_elements=len(prov_rows), total_bp=total_bp,
               pct_viral_cohort=100 * vbp / total_bp if total_bp else 0,
               enrichment_all=enr_all, phispy_consensus=consensus, checkv=dict(cv))
summary[f"enrichment_contigN50_ge_{N50_MIN}"] = enr_sens
json.dump(summary, open(os.path.join(RESULTS, "summary.json"), "w"), indent=2)

print(f"\n=== VDGCI SUMMARY: {len(accs)} genomes, {len(prov_rows)} viral elements, "
      f"{summary['pct_viral_cohort']:.2f}% of {total_bp/1e6:.1f} Mb ===")
for lab, e in [("all genomes", enr_all), (f"contig N50 >= {N50_MIN//1000} kb", enr_sens)]:
    if "fold" in e:
        print(f"defense enrichment ({lab}, n={e['n_genomes']}): {e['K']}/{e['S']} systems in viral DNA "
              f"({e['pct_systems_in_viral']:.1f}% vs {e['pct_sequence_viral']:.2f}% of sequence) = "
              f"{e['fold']:.2f}x, p={e['binom_p']:.2e}")
if consensus["pct"] is not None:
    print(f"PhiSpy corroboration (on {consensus['genomes_with_phispy']} annotated genomes): "
          f"{consensus['corroborated']}/{consensus['genomad_proviruses_in_those']} = {consensus['pct']:.0f}%")
if cv: print("CheckV tiers:", dict(cv))
print("Tables ->", RESULTS)
