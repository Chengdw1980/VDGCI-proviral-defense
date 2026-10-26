#!/bin/bash
# Build the three conda environments (once; ~30-60 min on a login node, no SUs)
# and write wrapper commands to $TOOLS_BIN so the pipeline finds the pinned
# DIAMOND 2.0.15 and DefenseFinder 3.0.0 ahead of anything else on PATH.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
cd "$PIPELINE_ROOT"
module load anaconda 2>/dev/null || true
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${PROJECT_BASE}/conda_pkgs}"
for pair in "env:environment.yml" "env_diamond:environment_diamond.yml" "env_df:environment_df.yml"; do
  name=${pair%%:*}; yml=${pair#*:}; prefix="${PROJECT_BASE}/vdgci/${name}"
  if [ -x "${prefix}/bin/python" ] || [ -x "${prefix}/bin/diamond" ]; then echo "[ok] ${name} exists"; continue; fi
  echo "[build] ${name} from ${yml} ..."
  conda env create -p "$prefix" -f "$yml"
done
mkdir -p "$TOOLS_BIN"
cat > "$TOOLS_BIN/diamond" <<W
#!/bin/bash
exec "${PROJECT_BASE}/vdgci/env_diamond/bin/diamond" "\$@"
W
cat > "$TOOLS_BIN/defense-finder" <<W
#!/bin/bash
export PATH="${PROJECT_BASE}/vdgci/env_df/bin:\$PATH"
exec "${PROJECT_BASE}/vdgci/env_df/bin/defense-finder" "\$@"
W
chmod +x "$TOOLS_BIN/diamond" "$TOOLS_BIN/defense-finder"
echo "Checks:"; "$TOOLS_BIN/diamond" version; "$TOOLS_BIN/defense-finder" run --help | head -1
echo "ALL ENVIRONMENTS READY"
