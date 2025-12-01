#!/bin/bash

# --------------------------------------------------------------------
# -- Optimized protein family clustering using MCL
# -- parallel post-processing
# --------------------------------------------------------------------

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Default parameters
INPUT_FILE='output/similarity_edgelists/filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv'
OUTPUT_DIR='output/clusters'
MCL_INFLATION=2.0
PRUNING_THRESHOLD=4000
SELECT_DOWN_TO=500
RECOVER=600
PCT=90
NUM_THREADS=$(nproc)
MIN_FAMILY_SIZE=2

# Parse arguments
while getopts "i:o:I:P:S:R:t:m:h" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        I) MCL_INFLATION="${OPTARG}" ;;
        P) PRUNING_THRESHOLD="${OPTARG}" ;;
        S) SELECT_DOWN_TO="${OPTARG}" ;;
        R) RECOVER="${OPTARG}" ;;
        t) NUM_THREADS="${OPTARG}" ;;
        m) MIN_FAMILY_SIZE="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE    Input edgelist file
  -o DIR     Output directory
  -I NUM     MCL inflation (default: 2.0, higher=more clusters)
  -P NUM     Pruning threshold (default: 4000)
  -S NUM     Select down to N entries (default: 500)
  -R NUM     Recover to N entries (default: 600)
  -t NUM     Number of threads for MCL (default: all CPUs)
  -m NUM     Minimum family size (default: 2)
  -h         Show this help

INFLATION GUIDELINES FOR Ks:
  1.4-2.0: Balanced families (default)
  2.5-3.0: Stricter clustering (breaks large families)
  3.5-5.0: Very strict (many small families)

EXAMPLES:
  # Standard clustering
  $0 -i edgelist.tsv -o clusters/

  # Strict clustering for Ks analysis
  $0 -i edgelist.tsv -I 2.5 -o clusters/

  # Very strict to break mega-families
  $0 -i edgelist.tsv -I 3.0 -P 2000 -S 300

EOF
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# Setup logging
LOG_DIR="logs/pipeline"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh)_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -i "$LOG_FILE") 2>&1

echo -e "${GREEN}===================================="
echo " OPTIMIZED MCL CLUSTERING"
echo "====================================${NC}"
echo " Input:     $INPUT_FILE"
echo " Output:    $OUTPUT_DIR"
echo " Inflation: $MCL_INFLATION"
echo " Threads:   $NUM_THREADS"
echo -e "${GREEN}====================================${NC}"

# Validate input
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}ERROR: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

# Check MCL installation
if ! command -v mcl &>/dev/null; then
    echo -e "${RED}ERROR: MCL not found. Install with:${NC}"
    echo "  conda install -c bioconda mcl"
    echo "  or apt-get install mcl"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build output filenames
BASE_NAME=$(basename "${INPUT_FILE%.tsv}")
OUTPUT_FILE="${OUTPUT_DIR}/protein_families_${BASE_NAME}.txt"
TSV_OUTPUT_FILE="${OUTPUT_FILE%.txt}.tsv"

# Count input edges
echo -e "${YELLOW}Analyzing input network...${NC}"
EDGE_COUNT=$(wc -l < "$INPUT_FILE")
NODE_COUNT=$(awk '{print $1; print $2}' "$INPUT_FILE" | sort -u | wc -l)
echo -e "${GREEN}✓ Nodes: $NODE_COUNT${NC}"
echo -e "${GREEN}✓ Edges: $EDGE_COUNT${NC}"

START_TIME=$(date +%s)

# ========================================
# STEP 1: Run MCL with optimizations
# ========================================
echo -e "${YELLOW}Running MCL clustering...${NC}"

# MCL with optimized parameters and threading
mcl "$INPUT_FILE" \
    -I "$MCL_INFLATION" \
    --abc \
    -V all \
    --discard-loops=y \
    -te "$NUM_THREADS" \
    -P "$PRUNING_THRESHOLD" \
    -S "$SELECT_DOWN_TO" \
    -R "$RECOVER" \
    -pct "$PCT" \
    -o "$OUTPUT_FILE" 2>&1 | while read line; do
        echo "  MCL: $line"
    done

MCL_END=$(date +%s)
MCL_TIME=$((MCL_END - START_TIME))
echo -e "${GREEN}✓ MCL completed in ${MCL_TIME}s${NC}"

# ========================================
# STEP 2: Post-process and fix merged IDs
# ========================================
echo -e "${YELLOW}Post-processing clusters...${NC}"

# Create optimized Python script for processing
cat > /tmp/process_clusters_$$.py <<'PYTHON_EOF'
#!/usr/bin/env python3
import sys
import re
from collections import defaultdict

input_file = sys.argv[1]
tsv_file = sys.argv[2]
min_size = int(sys.argv[3])

# Pattern to fix merged IDs (e.g., "ABC123DEF456" -> "ABC123 DEF456")
id_pattern = re.compile(r'([A-Z]+\d+)(?=[A-Z])')

# Read and process clusters
clusters = []
total_genes = 0
singleton_count = 0

with open(input_file) as f:
    for line in f:
        # Fix merged IDs
        fixed_line = id_pattern.sub(r'\1 ', line.strip())
        genes = fixed_line.split()
        
        if len(genes) >= min_size:
            clusters.append(genes)
            total_genes += len(genes)
        elif len(genes) == 1:
            singleton_count += 1

# Sort clusters by size (largest first)
clusters.sort(key=len, reverse=True)

# Write cleaned MCL format
with open(input_file, 'w') as f:
    for cluster in clusters:
        f.write(' '.join(cluster) + '\n')

# Write TSV format
with open(tsv_file, 'w') as f:
    f.write("geneName\tfamily\n")
    for family_id, cluster in enumerate(clusters, 1):
        for gene in cluster:
            f.write(f"{gene}\t{family_id}\n")

# Print statistics
print(f"Total clusters: {len(clusters)}")
print(f"Total genes: {total_genes}")
print(f"Singletons removed: {singleton_count}")
if clusters:
    sizes = [len(c) for c in clusters]
    print(f"Largest family: {max(sizes)} genes")
    print(f"Average family size: {sum(sizes)/len(sizes):.1f}")
    
    # Size distribution
    size_dist = defaultdict(int)
    for s in sizes:
        if s <= 10:
            size_dist[f"{s}"] += 1
        elif s <= 50:
            size_dist["11-50"] += 1
        elif s <= 100:
            size_dist["51-100"] += 1
        elif s <= 500:
            size_dist["101-500"] += 1
        else:
            size_dist["500+"] += 1
    
    print("\nFamily size distribution:")
    for size_range in ["2", "3", "4", "5", "6-10", "11-50", "51-100", "101-500", "500+"]:
        if size_range in size_dist or "-" in size_range:
            count = size_dist.get(size_range, 0)
            if "-" not in size_range:
                # Individual sizes 2-5
                count = size_dist.get(size_range, 0)
            if count > 0:
                print(f"  Size {size_range}: {count} families")
PYTHON_EOF

python3 /tmp/process_clusters_$$.py "$OUTPUT_FILE" "$TSV_OUTPUT_FILE" "$MIN_FAMILY_SIZE"

# Clean up
rm -f /tmp/process_clusters_$$.py

# ========================================
# STEP 3: Analyze results
# ========================================
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo ""
echo -e "${GREEN}===================================="
echo " CLUSTERING COMPLETE"
echo "====================================${NC}"
echo " Total time:     ${TOTAL_TIME}s"
echo " MCL time:       ${MCL_TIME}s"
echo " Post-process:   $((TOTAL_TIME - MCL_TIME))s"
echo -e "${GREEN}====================================${NC}"

# Final statistics
if [ -f "$TSV_OUTPUT_FILE" ]; then
    FAMILY_COUNT=$(tail -n +2 "$TSV_OUTPUT_FILE" | cut -f2 | sort -u | wc -l)
    CLUSTERED_GENES=$(tail -n +2 "$TSV_OUTPUT_FILE" | wc -l)
    
    echo ""
    echo "Clustering summary:"
    echo "  Total families: $FAMILY_COUNT"
    echo "  Clustered genes: $CLUSTERED_GENES"
    echo "  Average genes/family: $(echo "scale=1; $CLUSTERED_GENES / $FAMILY_COUNT" | bc)"
    
    # Check for problematic large families
    LARGEST=$(tail -n +2 "$TSV_OUTPUT_FILE" | cut -f2 | sort | uniq -c | sort -rn | head -1 | awk '{print $1}')
    
    if [ "$LARGEST" -gt 500 ]; then
        echo ""
        echo -e "${YELLOW}⚠ WARNING: Very large family detected (${LARGEST} genes)${NC}"
        echo "  For Ks analysis, consider:"
        echo "    1. Re-run with higher inflation: -I 2.5 or -I 3.0"
        echo "    2. Stricter BLAST filtering before clustering"
        echo "    3. Lower pruning threshold: -P 2000"
    fi
    
    # Estimate Ks computation load
    echo ""
    echo "Ks computation estimate:"
    TOTAL_PAIRS=$(tail -n +2 "$TSV_OUTPUT_FILE" | cut -f2 | sort | uniq -c | \
                  awk '{sum += $1*($1-1)/2} END {printf "%.0f", sum}')
    echo "  Total pairwise comparisons: $TOTAL_PAIRS"
    
    if [ "$TOTAL_PAIRS" -gt 10000000 ]; then
        echo -e "${YELLOW}  Note: This will generate >10M pairs for Ks calculation${NC}"
        echo "  Consider splitting large families or using stricter clustering"
    fi
fi

echo ""
echo "Output files:"
echo "  MCL format: $OUTPUT_FILE"
echo "  TSV format: $TSV_OUTPUT_FILE"
echo ""

exit 0