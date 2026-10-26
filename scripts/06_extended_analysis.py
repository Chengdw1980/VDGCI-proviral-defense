#!/usr/bin/env python3
"""
06_extended_analysis.py — pre-specified analyses for the 5,054-genome validation cohort
(VDGCI pipeline v3, October 2026).  Run AFTER 04_gtdb_tree.py and 05_analyze.py.

Reads the per-genome files kept in $KEEP/per_genome (same inputs as 05) plus the
pruned GTDB tree, and writes to $RESULTS/extended/:

  A. defense_systems_located.tsv   every defense system with its provirus (if any),
                                   that provirus's CheckV completeness/quality, and a
                                   mechanism class (see config/defense_mechanism_map.tsv)
  B. enrichment_by_mechanism.csv   fold-enrichment per mechanism class, genome-level
                                   bootstrap 95% CI, for three genome sets
                                   (all; contig N50 >= 20 kb; finished assemblies)
  C. enrichment_by_type.csv        per defense type: fold, binomial P, BH-FDR q
  D. enrichment_by_phylum.csv      per phylum with >= MIN_PHYLUM genomes, bootstrap CI
  E. decay_vs_defense.csv          defense systems per 100 kb of provirus by CheckV tier
                                   (tests whether defense cargo persists as prophages decay)
  F. phylo_signal.csv              Blomberg's K (permutation) + Pagel's lambda (ML/LRT)
                                   on the GTDB bac120 tree for per-genome traits
  G. genome_model.csv              log1p(% viral-derived) ~ log genome size + reduced-genome
                                   flag + fragmentation, phylum-clustered robust SEs
  H. mge_comparison.csv            defense enrichment in proviruses vs geNomad plasmids

Usage:  source config/config.sh && python3 scripts/06_extended_analysis.py
Optional env: N_BOOT (default 1000), MIN_PHYLUM (30), N_PERM (999), TREE (path)
"""
import os, sys, csv, glob, re, json
from collections import defaultdict
import numpy as np, warnings
warnings.filterwarnings("ignore")

def env(k, d=None):
    v = os.environ.get(k, d)
    if v is None: sys.exit(f"Missing {k}; run `source config/config.sh` first.")
    return v
KEEP = env("KEEP"); RESULTS = env("RESULTS"); ROOT = env("PIPELINE_ROOT")
PG = os.path.join(KEEP, "per_genome"); OUT = os.path.join(RESULTS, "extended"); os.makedirs(OUT, exist_ok=True)
N_BOOT = int(os.environ.get("N_BOOT", "1000")); MIN_PHYLUM = int(os.environ.get("MIN_PHYLUM", "30"))
N_PERM = int(os.environ.get("N_PERM", "999")); N50_MIN = int(os.environ.get("SENS_MIN_N50", "20000"))
TREE = os.environ.get("TREE", os.path.join(os.environ.get("PHYLO_OUT", ""), "cohort_gtdb_bac120.nwk"))
rng = np.random.default_rng(int(os.environ.get("SEED", "20261023")))

# ---------------- metadata ----------------
meta = {}
for cand in [os.environ.get("MANIFEST", ""), os.path.join(ROOT, "config", "cohort_manifest.tsv")]:
    if cand and os.path.exists(cand):
        for r in csv.DictReader(open(cand), delimiter="\t"): meta.setdefault(r["accession"], r)
asm = {}
p = os.path.join(ROOT, "config", "assembly_stats_HQ.tsv")
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter="\t"): asm[r["accession"]] = r
mech_map = {}
p = os.path.join(ROOT, "config", "defense_mechanism_map.tsv")
if os.path.exists(p):
    for r in csv.DictReader(open(p), delimiter="\t"): mech_map[r["type"]] = r["mechanism_class"]
def mech(t):
    if t in mech_map: return mech_map[t]
    if t.startswith("RM_"): return "Direct nucleic-acid targeting"
    if t.startswith("Abi"): return "Abortive infection / toxin-antitoxin"
    return "Other or mechanism unresolved"
# Pre-specified genome-reduced / host-restricted lineages (GTDB r232 names):
REDUCED_PHYLA = {"Patescibacteriota", "Chlamydiota"}
REDUCED_ORDERS = {"Mycoplasmatales", "Rickettsiales", "Holosporales"}
def is_reduced(g): return g["phylum"] in REDUCED_PHYLA or g.get("order", "") in REDUCED_ORDERS

# ---------------- per-genome parsing ----------------
G = {}          # acc -> dict
PROV = defaultdict(list)   # (acc, contig) -> [(s, e, pid)]
PLAS = defaultdict(list)   # (acc, contig) -> [(1, L)] whole plasmid contigs
SYS = []        # rows
CV = {}         # (acc, contig, s, e) -> (completeness, quality)

def parse_cv_id(cid, L=None):
    contig = cid.split("|")[0]
    m = re.search(r"provirus_(\d+)_(\d+)", cid)
    if m: return contig, int(m.group(1)), int(m.group(2))
    return contig, 1, int(L) if L else None

for d in sorted(glob.glob(os.path.join(PG, "*"))):
    acc = os.path.basename(d); vs = os.path.join(d, f"{acc}_virus_summary.tsv")
    if not os.path.exists(vs): continue
    gs = {"genome_bp": 0, "n_contigs": 0, "contig_n50": 0}
    gsf = os.path.join(d, "genome_stats.tsv")
    if os.path.exists(gsf):
        r = next(csv.DictReader(open(gsf), delimiter="\t")); gs = {k: int(r[k]) for k in gs}
    m = meta.get(acc, {})
    g = dict(acc=acc, phylum=m.get("phylum", ""), order=m.get("order", ""), family=m.get("family", ""), viral_bp=0, n_prov=0,
             plasmid_bp=0, level=(asm.get(acc, {}).get("assembly_level") or m.get("ncbi_assembly_level", "")), **gs)
    for r in csv.DictReader(open(vs), delimiter="\t"):
        sn = r["seq_name"]; L = int(r["length"]); co = (r.get("coordinates") or "").strip()
        if co and co.upper() != "NA" and "-" in co:
            contig = sn.split("|")[0]; s, e = map(int, co.split("-"))
        else: contig, s, e = sn, 1, L
        PROV[(acc, contig)].append((s, e, f"{contig}:{s}-{e}")); g["viral_bp"] += e - s + 1; g["n_prov"] += 1
    ps = os.path.join(d, f"{acc}_plasmid_summary.tsv")
    if os.path.exists(ps):
        for r in csv.DictReader(open(ps), delimiter="\t"):
            L = int(r["length"]); PLAS[(acc, r["seq_name"])].append((1, L)); g["plasmid_bp"] += L
    q = os.path.join(d, "checkv_quality_summary.tsv")
    if os.path.exists(q):
        for r in csv.DictReader(open(q), delimiter="\t"):
            c, s, e = parse_cv_id(r["contig_id"], r.get("contig_length"))
            try: comp = float(r.get("completeness", "nan"))
            except ValueError: comp = float("nan")
            CV[(acc, c, s, e)] = (comp, r.get("checkv_quality", "NA"))
    gene_xy = {}
    hf = os.path.join(d, f"{acc}_protein_headers.txt")
    if os.path.exists(hf):
        for line in open(hf):
            pp = line[1:].split("#")
            try: gene_xy[pp[0].strip()] = (pp[0].strip().rsplit("_", 1)[0], int(pp[1]), int(pp[2]))
            except (IndexError, ValueError): pass
    sf = glob.glob(os.path.join(d, "*defense_finder_systems.tsv")); gf = glob.glob(os.path.join(d, "*defense_finder_genes.tsv"))
    g["has_df"] = bool(sf and gf)
    if g["has_df"]:
        sg = defaultdict(list)
        for r in csv.DictReader(open(gf[0]), delimiter="\t"): sg[r.get("sys_id", "")].append(r.get("hit_id", ""))
        sg = {k: list(dict.fromkeys(v)) for k, v in sg.items()}; seen = set()  # de-duplicate (v3.2)
        for r in csv.DictReader(open(sf[0]), delimiter="\t"):
            if r.get("sys_id", "") in seen: continue  # DefenseFinder can list a system twice
            seen.add(r.get("sys_id", ""))
            pid_hit, pl_hit = None, False
            for hid in sg.get(r.get("sys_id", ""), []):
                xy = gene_xy.get(hid)
                if not xy: continue
                c, s, e = xy
                for ps_, pe_, pid in PROV.get((acc, c), []):
                    if s <= pe_ and e >= ps_: pid_hit = pid
                if PLAS.get((acc, c)): pl_hit = True
            comp, qual = (float("nan"), "")
            if pid_hit:
                c, rng_ = pid_hit.split(":"); s, e = map(int, rng_.split("-"))
                comp, qual = CV.get((acc, c, s, e), (float("nan"), "NA"))
            t = r.get("type", "")
            SYS.append(dict(acc=acc, sys_id=r.get("sys_id", ""), type=t, subtype=r.get("subtype", ""),
                            mechanism=mech(t), provirus=pid_hit or "", in_viral=bool(pid_hit),
                            in_plasmid=pl_hit and not pid_hit, checkv_completeness=comp, checkv_quality=qual))
    G[acc] = g
if not G: sys.exit(f"No per-genome results in {PG}")
accs = sorted(a for a in G if G[a]["has_df"])
print(f"{len(G)} genomes parsed; {len(accs)} with DefenseFinder; {len(SYS)} systems")

def write(name, rows, header):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t" if name.endswith(".tsv") else ","); w.writerow(header); w.writerows(rows)

# A
write("defense_systems_located.tsv", [[r[k] for k in ("acc", "sys_id", "type", "subtype", "mechanism", "provirus",
      "in_viral", "in_plasmid", "checkv_completeness", "checkv_quality")] for r in SYS],
      ["accession", "sys_id", "type", "subtype", "mechanism_class", "provirus", "in_viral_dna", "in_plasmid",
       "provirus_checkv_completeness", "provirus_checkv_quality"])

# ---------------- genome-clustered enrichment ----------------
def per_genome_counts(filter_fn, which="in_viral"):
    k = defaultdict(int); n = defaultdict(int)
    for r in SYS:
        if filter_fn(r): n[r["acc"]] += 1; k[r["acc"]] += r[which]
    return k, n
def fold_ci(sub, filter_fn=lambda r: True, which="in_viral", bp_key="viral_bp"):
    sub = [a for a in sub if a in G]
    if not sub: return None
    k, n = per_genome_counts(filter_fn, which)
    K = np.array([k[a] for a in sub], float); N = np.array([n[a] for a in sub], float)
    V = np.array([G[a][bp_key] for a in sub], float); B = np.array([G[a]["genome_bp"] for a in sub], float)
    def f(ix):
        if N[ix].sum() == 0 or V[ix].sum() == 0: return np.nan
        return (K[ix].sum() / N[ix].sum()) / (V[ix].sum() / B[ix].sum())
    est = f(np.arange(len(sub)))
    bs = np.array([f(rng.integers(0, len(sub), len(sub))) for _ in range(N_BOOT)]); bs = bs[np.isfinite(bs)]
    lo, hi = (np.percentile(bs, [2.5, 97.5]) if len(bs) else (np.nan, np.nan))
    return dict(n_genomes=len(sub), n_systems=int(N.sum()), n_in=int(K.sum()),
                pct_seq=100 * V.sum() / B.sum() if B.sum() else np.nan, fold=est, lo=lo, hi=hi)

sets = {"all": accs, f"contigN50>={N50_MIN}": [a for a in accs if G[a]["contig_n50"] >= N50_MIN],
        "finished_assemblies": [a for a in accs if G[a]["level"] in ("Complete Genome", "Chromosome")]}
MECHS = sorted({r["mechanism"] for r in SYS})
rows = []
for sname, sub in sets.items():
    for mname in ["ALL"] + MECHS:
        ff = (lambda r: True) if mname == "ALL" else (lambda r, m=mname: r["mechanism"] == m)
        e = fold_ci(sub, ff)
        if e: rows.append([sname, mname, e["n_genomes"], e["n_systems"], e["n_in"], round(e["pct_seq"], 3),
                           round(e["fold"], 3), round(e["lo"], 3), round(e["hi"], 3)])
write("enrichment_by_mechanism.csv", rows, ["genome_set", "mechanism_class", "n_genomes", "n_systems",
      "n_in_viral_dna", "pct_sequence_viral", "fold", "boot95_lo", "boot95_hi"])
# Pre-specified contrast: Abi/TA vs direct, as a ratio of folds with a bootstrap CI
def contrast(sub):
    k1, n1 = per_genome_counts(lambda r: r["mechanism"].startswith("Abortive"))
    k2, n2 = per_genome_counts(lambda r: r["mechanism"].startswith("Direct"))
    A = np.array([[k1[a], n1[a], k2[a], n2[a]] for a in sub], float)
    def rr(ix):
        s = A[ix].sum(0)
        return (s[0] / s[1]) / (s[2] / s[3]) if s[1] and s[3] and s[2] else np.nan
    est = rr(np.arange(len(sub))); bs = np.array([rr(rng.integers(0, len(sub), len(sub))) for _ in range(N_BOOT)])
    bs = bs[np.isfinite(bs)]
    return est, *(np.percentile(bs, [2.5, 97.5]) if len(bs) else (np.nan, np.nan))
con = {s: [round(x, 3) for x in contrast(sub)] for s, sub in sets.items()}

# C. per type with BH-FDR (binomial; descriptive — clustering handled in B/D)
from scipy import stats
vb = sum(G[a]["viral_bp"] for a in accs); gb = sum(G[a]["genome_bp"] for a in accs); p0 = vb / gb if gb else 0
bt = defaultdict(lambda: [0, 0])
for r in SYS:
    if r["acc"] in G: bt[r["type"]][0] += 1; bt[r["type"]][1] += r["in_viral"]
tr = []
for t, (n, k) in bt.items():
    tr.append([t, mech(t), n, k, round(100 * k / n, 2), round((k / n) / p0, 3) if p0 else "",
               stats.binomtest(k, n, p0, alternative="greater").pvalue if p0 else 1.0])
if tr:
    pv = np.array([x[-1] for x in tr]); o = np.argsort(pv); m = len(pv)
    q = np.empty(m); q[o] = np.minimum.accumulate((pv[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    for x, qq in zip(tr, q): x.append(min(qq, 1.0))
write("enrichment_by_type.csv", sorted(tr, key=lambda x: -x[2]),
      ["type", "mechanism_class", "n_systems", "n_in_viral_dna", "pct_in_viral_dna", "fold", "binom_p", "bh_q"])

# D. per phylum
byp = defaultdict(list)
for a in accs: byp[G[a]["phylum"]].append(a)
rows = []
for ph, sub in sorted(byp.items(), key=lambda kv: -len(kv[1])):
    if len(sub) < MIN_PHYLUM: continue
    e = fold_ci(sub)
    if e: rows.append([ph, e["n_genomes"], e["n_systems"], e["n_in"], round(e["pct_seq"], 3),
                       round(e["fold"], 3) if np.isfinite(e["fold"]) else "", round(e["lo"], 3), round(e["hi"], 3)])
write("enrichment_by_phylum.csv", rows, ["phylum", "n_genomes", "n_systems", "n_in_viral_dna",
      "pct_sequence_viral", "fold", "boot95_lo", "boot95_hi"])

# E. decay: systems per 100 kb of provirus, by CheckV tier
tier_len = defaultdict(float); tier_sys = defaultdict(int); tier_of = {}
for (a, c), ivs in PROV.items():
    if a not in G or not G[a]["has_df"]: continue
    for s, e, pid in ivs:
        comp, qual = CV.get((a, c, s, e), (float("nan"), "NA"))
        tier = qual if qual in ("Complete", "High-quality", "Medium-quality", "Low-quality") else "Not-determined"
        if tier == "Complete": tier = "High-quality"
        tier_len[tier] += e - s + 1; tier_of[(a, pid)] = tier
for r in SYS:
    if r["in_viral"]: tier_sys[tier_of.get((r["acc"], r["provirus"]), "Not-determined")] += 1
rows = [[t, tier_sys[t], round(tier_len[t]), round(1e5 * tier_sys[t] / tier_len[t], 4) if tier_len[t] else ""]
        for t in ("High-quality", "Medium-quality", "Low-quality", "Not-determined")]
write("decay_vs_defense.csv", rows, ["checkv_tier", "n_defense_systems", "provirus_bp", "systems_per_100kb"])

# H. MGE comparison (plasmids identified by geNomad)
mge = []
for lab, which, key in [("proviruses", "in_viral", "viral_bp"), ("geNomad plasmids", "in_plasmid", "plasmid_bp")]:
    e = fold_ci(accs, which=which, bp_key=key)
    if e: mge.append([lab, e["n_genomes"], e["n_systems"], e["n_in"], round(e["pct_seq"], 3),
                      round(e["fold"], 3) if np.isfinite(e["fold"]) else "", round(e["lo"], 3), round(e["hi"], 3)])
write("mge_comparison.csv", mge, ["compartment", "n_genomes", "n_systems", "n_in_compartment",
      "pct_sequence", "fold", "boot95_lo", "boot95_hi"])

# ---------------- F. phylogenetic signal ----------------
def vcv_fast(tree, tips):
    import treeswift
    t = tree.extract_tree_with(set(tips)); t.suppress_unifurcations()
    idx = {lab: i for i, lab in enumerate(tips)}; n = len(tips); C = np.zeros((n, n)); depth = {}
    for nd in t.traverse_preorder():
        depth[nd] = 0.0 if nd.is_root() else depth[nd.parent] + (nd.edge_length or 0.0)
    leaves = {}
    for nd in t.traverse_postorder():
        if nd.is_leaf(): leaves[nd] = [idx[nd.label]]; C[idx[nd.label], idx[nd.label]] = depth[nd]; continue
        ch = [leaves[c] for c in nd.children]
        for i in range(len(ch)):
            for j in range(i + 1, len(ch)):
                C[np.ix_(ch[i], ch[j])] = depth[nd]; C[np.ix_(ch[j], ch[i])] = depth[nd]
        leaves[nd] = [x for c in ch for x in c]
        for c in nd.children: del leaves[c]
    return C / C.diagonal().max()
def blomberg_k(x, C, nperm):
    n = len(x); Ci = np.linalg.inv(C); one = np.ones(n); s1 = one @ Ci @ one
    exp = (np.trace(C) - n / s1) / (n - 1)
    def K(xx):
        a = (one @ Ci @ xx) / s1; r = xx - a
        return ((r @ r) / (r @ Ci @ r)) / exp
    k = K(x); null = np.array([K(rng.permutation(x)) for _ in range(nperm)])
    return k, (1 + (null >= k).sum()) / (nperm + 1)
def pagel_lambda(x, C):
    from scipy import optimize
    n = len(x); D = np.diag(np.diag(C)); one = np.ones(n)
    def nll(lam):
        V = lam * C + (1 - lam) * D; sign, ld = np.linalg.slogdet(V); Vi = np.linalg.inv(V)
        a = (one @ Vi @ x) / (one @ Vi @ one); r = x - a; s2 = (r @ Vi @ r) / n
        return 0.5 * (n * np.log(2 * np.pi * s2) + ld + n)
    res = optimize.minimize_scalar(nll, bounds=(0, 1), method="bounded", options={"xatol": 1e-3})
    lr = 2 * (nll(0.0) - res.fun)
    return res.x, lr, 0.5 * stats.chi2.sf(max(lr, 0), 1)
phylo_rows = []
if os.path.exists(TREE):
    import treeswift
    tree = treeswift.read_tree_newick(TREE)
    tips = [l.label for l in tree.traverse_leaves() if l.label in G]
    if len(tips) >= 20:
        C = vcv_fast(tree, tips)
        ndef = defaultdict(int)
        for r in SYS: ndef[r["acc"]] += 1
        traits = {"pct_viral_log1p": [np.log1p(100 * G[a]["viral_bp"] / G[a]["genome_bp"]) if G[a]["genome_bp"] else np.nan for a in tips],
                  "n_proviruses_log1p": [np.log1p(G[a]["n_prov"]) for a in tips],
                  "n_defense_log1p": [np.log1p(ndef[a]) if G[a]["has_df"] else np.nan for a in tips],
                  "genome_size_log (control)": [np.log(G[a]["genome_bp"]) if G[a]["genome_bp"] else np.nan for a in tips]}
        for name, x in traits.items():
            x = np.array(x, float); ok = np.isfinite(x)
            Cs = C[np.ix_(ok, ok)]; k, pk = blomberg_k(x[ok], Cs, N_PERM); lam, lr, pl = pagel_lambda(x[ok], Cs)
            phylo_rows.append([name, int(ok.sum()), round(k, 4), pk, round(lam, 4), round(lr, 3), pl])
            print(f"  phylo signal {name}: K={k:.3f} (P={pk:.3g}) lambda={lam:.3f} (P={pl:.3g})")
    import shutil; shutil.copy(TREE, os.path.join(OUT, "cohort_gtdb_bac120.nwk"))
else:
    print(f"  (no tree at {TREE}; run scripts/04_gtdb_tree.py first — phylogenetic signal skipped)")
write("phylo_signal.csv", phylo_rows, ["trait", "n", "blomberg_K", "P_K_perm", "pagel_lambda", "LR", "P_lambda"])

# ---------------- G. genome-level model ----------------
try:
    import statsmodels.api as sm
    X, y, grp = [], [], []
    for a in G:
        g = G[a]
        if not g["genome_bp"]: continue
        y.append(np.log1p(100 * g["viral_bp"] / g["genome_bp"]))
        X.append([np.log(g["genome_bp"] / 1e6), float(is_reduced(g)),
                  float(g["contig_n50"] < N50_MIN), float(g["level"] in ("Complete Genome", "Chromosome"))])
        grp.append(g["phylum"] or "NA")
    X = sm.add_constant(np.array(X)); codes = {p: i for i, p in enumerate(sorted(set(grp)))}
    fit = sm.OLS(np.array(y), X).fit(cov_type="cluster", cov_kwds={"groups": np.array([codes[p] for p in grp])})
    names = ["intercept", "log_genome_Mb", "reduced_lineage", "fragmented_N50<20kb", "finished_assembly"]
    write("genome_model.csv", [[n_, round(b, 4), round(se, 4), p_] for n_, b, se, p_ in
          zip(names, fit.params, fit.bse, fit.pvalues)], ["term", "coef", "cluster_robust_se", "P"])
except Exception as ex:
    print("  genome model skipped:", ex)

summary = dict(n_genomes=len(G), n_with_defensefinder=len(accs), n_systems=len(SYS),
               abiTA_vs_direct_fold_ratio={k: dict(est=v[0], lo=v[1], hi=v[2]) for k, v in con.items()})
json.dump(summary, open(os.path.join(OUT, "extended_summary.json"), "w"), indent=2, default=float)
print("Abi/TA vs direct (ratio of folds, est/lo/hi):", con)
print("Extended tables ->", OUT)
