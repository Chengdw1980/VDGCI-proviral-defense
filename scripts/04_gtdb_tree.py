#!/usr/bin/env python3
"""
Host phylogeny for the cohort, taken from the published GTDB r232 bacterial
reference tree (bac120: 120 concatenated single-copy marker proteins), pruned
to the genomes in $ACCESSION_LIST.

Why not a 16S tree: most draft genomes / MAGs in the GTDB family cohort have no
assembled 16S gene, so a 16S tree would silently drop them. Every cohort genome
is a GTDB species representative and therefore a tip of the bac120 tree.

Tip labels in GTDB are 'RS_GCF_xxx' / 'GB_GCA_xxx'; we strip the prefix.
Output: $PHYLO_OUT/cohort_gtdb_bac120.nwk  (+ a list of any accessions not found)
Usage:  source config/config.sh && python3 scripts/04_gtdb_tree.py [tree] [accessions]
"""
import os, sys, warnings, treeswift
warnings.filterwarnings("ignore")

tree_path = sys.argv[1] if len(sys.argv) > 1 else os.environ["GTDB_TREE"]
acc_path  = sys.argv[2] if len(sys.argv) > 2 else os.environ["ACCESSION_LIST"]
out_dir   = os.environ.get("PHYLO_OUT", ".")
os.makedirs(out_dir, exist_ok=True)

want = {l.strip() for l in open(acc_path) if l.strip()}
t = treeswift.read_tree_newick(tree_path)
label = {}
for leaf in t.traverse_leaves():
    lab = leaf.label or ""
    acc = lab[3:] if lab[:3] in ("RS_", "GB_") else lab
    label[acc] = lab
present = {label[a] for a in want if a in label}
missing = sorted(a for a in want if a not in label)
pruned = t.extract_tree_with(present)
pruned.suppress_unifurcations()
for leaf in pruned.traverse_leaves():
    if leaf.label[:3] in ("RS_", "GB_"):
        leaf.label = leaf.label[3:]
out = os.path.join(out_dir, "cohort_gtdb_bac120.nwk")
pruned.write_tree_newick(out)
with open(os.path.join(out_dir, "cohort_gtdb_missing.txt"), "w") as f:
    f.write("\n".join(missing) + ("\n" if missing else ""))
print(f"GTDB tree: {len(present)}/{len(want)} cohort genomes placed -> {out}; missing {len(missing)}")
