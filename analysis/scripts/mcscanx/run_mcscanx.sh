#!/usr/bin/env bash
set -euo pipefail

# Run MCScanX for synteny/WGD detection using prepared inputs
# Usage:
#   bash analysis/run_mcscanx.sh output/mcscanx/soybean

if [[ ${1:-} == "" ]]; then
  echo "Usage: $0 <prefix_path_without_extension>"
  echo "  Example: $0 output/mcscanx/soybean"
  exit 1
fi

PREFIX="$1"
DIR=$(dirname "$PREFIX")
BASE=$(basename "$PREFIX")

GFF="${PREFIX}.gff"
BLAST="${PREFIX}.blast"

# Check inputs
if [[ ! -f "$GFF" ]]; then
  echo "ERROR: GFF not found: $GFF" >&2
  exit 2
fi
if [[ ! -f "$BLAST" ]]; then
  echo "ERROR: BLAST not found: $BLAST" >&2
  exit 2
fi

# Check MCScanX installed
if ! command -v MCScanX >/dev/null 2>&1; then
  echo "MCScanX not found. Install via conda:"
  echo "  conda install -c bioconda mcscanx"
  exit 3
fi

pushd "$DIR" >/dev/null
  echo "Running MCScanX on prefix: $BASE"
  MCScanX "$BASE"
  echo "\nMCScanX finished. Outputs:"
  ls -lh ${BASE}.* | sed 's/^/  - /'

  # Optional post-processing: classify duplicates if file exists
  if [[ -f "${BASE}.collinearity" ]]; then
    echo "\nParsing anchors and annotating with Ks (if available) ..."
    # Get absolute path to repo root
    REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    python3 "${REPO_ROOT}/scripts/annotate_mcscanx_anchors_with_ks.py" \
      --col ${BASE}.collinearity \
      --ks "${REPO_ROOT}/output/ks_results/ks_results_filtered.tsv" \
      --out ${BASE}_anchors_with_ks.tsv \
      --plots ${BASE}_anchors_plots
  fi
popd >/dev/null

echo "\n✓ MCScanX workflow complete."
