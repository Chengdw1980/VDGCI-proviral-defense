#!/bin/bash
# Build a host 16S rRNA phylogeny for the whole cohort (single job).
# Extracts one 16S gene per genome from its GBFF, aligns (MAFFT), infers ML tree
# (FastTree, GTR). Output: $PHYLO_OUT/cohort_16S.nwk
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"

python3 - <<'PY'
import os, glob
from Bio import SeqIO
GBFF_DIR=os.environ["GBFF_DIR"]; PHYLO_OUT=os.environ["PHYLO_OUT"]
out=os.path.join(PHYLO_OUT,"cohort_16S.fasta"); n=0
with open(out,"w") as w:
    for gb in sorted(glob.glob(os.path.join(GBFF_DIR,"*.gbff"))):
        acc=os.path.basename(gb)[:-5]; best=None
        for rec in SeqIO.parse(gb,"genbank"):
            for feat in rec.features:
                if feat.type=="rRNA":
                    prod=" ".join(feat.qualifiers.get("product",[""])).lower()
                    if "16s" in prod:
                        seq=str(feat.extract(rec.seq))
                        if best is None or len(seq)>len(best): best=seq
        if best:
            w.write(f">{acc}\n{best}\n"); n+=1
        else:
            print("no 16S in",acc)
print("wrote",n,"16S sequences ->",out)
PY

cd "$PHYLO_OUT"
mafft --auto --thread "$THREADS" cohort_16S.fasta > cohort_16S.aln 2> "${LOGS}/mafft.log"
FastTree -nt -gtr cohort_16S.aln > cohort_16S.nwk 2> "${LOGS}/fasttree.log"
echo "Tree written: ${PHYLO_OUT}/cohort_16S.nwk"
