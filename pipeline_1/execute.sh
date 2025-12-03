#!/bin/bash
# Master script to run all Pipeline 1 steps with progress bars and pretty output

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
NC='\033[0m'


# Usage/help
if [[ "$*" =~ -h ]]; then
  cat <<EOF
Usage: $0 [-s SPECIES]

OPTIONS:
  -s SPECIES   Input (species name, passed to all steps)
  -h           Show this help

All other arguments are forwarded to each step script.
EOF
  exit 0
fi

# Parse -s (and collect extra args)
SPECIES=""
EXTRA_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s) SPECIES="$2"; shift 2;;
    -h) shift;;
    *) EXTRA_ARGS+=("$1"); shift;;
  esac
done

# Autodetect BLAST results and set START_FROM automatically
START_FROM=1
if [ -n "$SPECIES" ]; then
  BLAST_DIR="$(cd "$(dirname "$0")/.." && pwd)/output/pipeline1/${SPECIES}/blast_results"
  if [ -f "$BLAST_DIR/blast_results_with_coverage.tsv" ]; then
    echo -e "${YELLOW}BLAST results with coverage found for species '$SPECIES'.${NC}"
    echo -e "=== ${YELLOW}Starting from step 4${NC} ===\n"
    START_FROM=4
  elif [ -f "$BLAST_DIR/blast_results.tsv" ]; then
    echo -e "${YELLOW}BLAST results found for species '$SPECIES'.${NC}"
    echo -e "=== ${YELLOW}Starting from step 3${NC} ===\n"
    START_FROM=3
  fi
fi

STEPS=(
  "1_filter_isoforms.sh:Filtering isoforms"
  "2_blast.sh:Running BLAST"
  "3_compute_coverage.sh:Computing coverage"
  "4_filter_pairs.sh:Filtering pairs"
  "5_prepare_edgelist.sh:Preparing edgelist"
  "6_cluster_families.sh:Clustering families"
)

function progress_bar() {
  local duration=${1:-5}
  local step_msg="${2:-Processing}"
  local bar_length=30
  echo -ne "${BLUE}${step_msg}... "
  for ((i=0; i<=bar_length; i++)); do
    percent=$((i * 100 / bar_length))
    bar=$(printf '%0.s#' $(seq 1 $i))
    spaces=$(printf '%0.s ' $(seq 1 $((bar_length-i))))
    echo -ne "[${bar}${spaces}] ${percent}%\r"
    sleep $(echo "scale=2; $duration/$bar_length" | bc)
  done
  echo -ne "${NC}"
}

echo -e "${GREEN}===================================="
echo " PIPELINE 1: Automated Execution"
echo -e "====================================${NC}"


step_num=1
for idx in "${!STEPS[@]}"; do
  step_idx=$((idx+1))
  if [ "$step_idx" -lt "$START_FROM" ]; then
    continue
  fi
  step="${STEPS[$idx]}"
  script="${step%%:*}"
  msg="${step#*:}"
  echo -e "${YELLOW}Step $step_idx: $msg${NC}"
  progress_bar 2 "$msg"
  # Build argument list for each step
  step_args=()
  [ -n "$SPECIES" ] && step_args+=("-s" "$SPECIES")
  step_args+=("${EXTRA_ARGS[@]}")
  if bash "$(dirname "$0")/$script" "${step_args[@]}"; then
    echo -e "${GREEN}✓ $msg completed${NC}\n"
  else
    echo -e "${RED}✗ $msg failed. Stopping pipeline.${NC}"
    exit 1
  fi
done

echo -e "${GREEN}===================================="
echo " PIPELINE 1: All steps finished! "
echo -e "====================================${NC}"
