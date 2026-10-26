#!/usr/bin/env python3
"""After the 10-genome TEST, measure real cost and project the full run.
Reads $KEEP/done/<ACC> files (accession, wall seconds, threads).
Usage: source config/config.sh && python3 scripts/07_cost_report.py [N_target]"""
import os, sys, glob, statistics as st
keep = os.environ["KEEP"]; n_target = int(sys.argv[1]) if len(sys.argv) > 1 else 5054
rows = []
for f in glob.glob(os.path.join(keep, "done", "*")):
    p = open(f).read().split()
    if len(p) >= 3: rows.append((p[0], int(p[1]), int(p[2])))
if not rows: sys.exit("No finished genomes yet.")
su = [w / 3600 * t for _, w, t in rows]          # cores x hours
mean, med, mx = st.mean(su), st.median(su), max(su)
proj = mean * n_target * 1.15                      # +15% for chunk tail idle time
print(f"finished genomes: {len(rows)}")
print(f"SUs per genome: mean {mean:.2f}, median {med:.2f}, max {mx:.2f}  (wall mean {st.mean(w for _,w,_ in rows)/60:.0f} min at {rows[0][2]} threads)")
print(f"projected for {n_target} genomes: ~{proj:,.0f} SUs (includes 15% packing overhead)")
print("Allocation on Anvil CPU: 50,000 SUs ->", "OK" if proj < 40000 else "TIGHT - review before full run")
