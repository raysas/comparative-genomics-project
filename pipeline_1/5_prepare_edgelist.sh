#!/bin/bash

# --------------------------------------------------------------------
# -- Optimized edgelist generation from BLAST results
# -- parallel processing and efficient sorting
# --------------------------------------------------------------------

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Default parameters
INPUT_FILE='output/blast_filtered/filtered_blast_results_id30_qcov50_scov50.tsv'
OUTPUT_DIR='output/similarity_edgelists'
WEIGHT_COLUMN_INDEX=12  # Bit score column
NUM_JOBS=$(nproc)
CHUNK_SIZE=10000000  # 10M lines per chunk
USE_MEMORY_SORT=false

# Parse arguments
while getopts "i:o:w:j:c:mh" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        w) WEIGHT_COLUMN_INDEX="${OPTARG}" ;;
        j) NUM_JOBS="${OPTARG}" ;;
        c) CHUNK_SIZE="${OPTARG}" ;;
        m) USE_MEMORY_SORT=true ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE    Input filtered BLAST results
  -o DIR     Output directory
  -w NUM     Weight column index (default: 12 for bit score)
  -j NUM     Parallel jobs (default: all CPUs)
  -c NUM     Chunk size for processing (default: 10000000)
  -m         Use memory-based sorting (faster for <100M edges)
  -h         Show this help

EXAMPLES:
  # Standard usage
  $0 -i blast_filtered.tsv -o edgelists/

  # Fast processing with memory sort
  $0 -i blast_filtered.tsv -m -j 32

  # Large file with chunking
  $0 -i huge_blast.tsv -c 50000000 -j 64

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
echo " OPTIMIZED EDGELIST GENERATION"
echo "====================================${NC}"
echo " Input:  $INPUT_FILE"
echo " Output: $OUTPUT_DIR"
echo " Weight column: $WEIGHT_COLUMN_INDEX"
echo " CPUs:   $NUM_JOBS"
if [ "$USE_MEMORY_SORT" = true ]; then
    echo " Mode:   Memory-based sorting"
else
    echo " Mode:   Disk-based sorting"
fi
echo -e "${GREEN}====================================${NC}"

# Validate inputs
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}ERROR: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

# Validate weight column is numeric
if ! [[ "$WEIGHT_COLUMN_INDEX" =~ ^[0-9]+$ ]]; then
    echo -e "${YELLOW}Warning: Weight column must be numeric. Using default (12)${NC}"
    WEIGHT_COLUMN_INDEX=12
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build output filename
OUTPUT_FILE="$OUTPUT_DIR/$(basename "${INPUT_FILE%.tsv}_wcol${WEIGHT_COLUMN_INDEX}_network.tsv")"

# Count input lines
echo -e "${YELLOW}Analyzing input file...${NC}"
TOTAL_LINES=$(wc -l < "$INPUT_FILE")
echo -e "${GREEN}✓ Input lines: $TOTAL_LINES${NC}"

# Create temp directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

START_TIME=$(date +%s)

# ========================================
# STEP 1: Extract and normalize edges
# ========================================
echo -e "${YELLOW}Step 1: Extracting edges and normalizing...${NC}"

if [ "$TOTAL_LINES" -lt 50000000 ]; then
    # Small file: single AWK pass
    echo "  Using single-pass processing..."
    
    awk -v wcol="$WEIGHT_COLUMN_INDEX" '
    NR > 1 {  # Skip header if present
        # Skip self-loops immediately
        if ($1 == $2) next
        
        # Normalize edge order (smaller ID first)
        if ($1 < $2) {
            print $1, $2, $(wcol)
        } else {
            print $2, $1, $(wcol)
        }
    }' "$INPUT_FILE" > "$TEMP_DIR/normalized.txt"
    
else
    # Large file: parallel chunk processing
    echo "  Using parallel chunk processing..."
    
    # Split into chunks
    tail -n +2 "$INPUT_FILE" | split -l "$CHUNK_SIZE" - "$TEMP_DIR/chunk_"
    
    # Function to process a chunk
    process_chunk() {
        local chunk="$1"
        local wcol="$2"
        local output="${chunk}.normalized"
        
        awk -v wcol="$wcol" '
        {
            if ($1 == $2) next
            if ($1 < $2) {
                print $1, $2, $(wcol)
            } else {
                print $2, $1, $(wcol)
            }
        }' "$chunk" > "$output"
        
        echo "$output"
    }
    
    export -f process_chunk
    
    # Process chunks in parallel
    ls "$TEMP_DIR"/chunk_* | \
        parallel -j "$NUM_JOBS" --progress process_chunk {} "$WEIGHT_COLUMN_INDEX" > "$TEMP_DIR/chunk_list.txt"
    
    # Combine chunks
    cat $(cat "$TEMP_DIR/chunk_list.txt") > "$TEMP_DIR/normalized.txt"
    
    # Clean up chunks
    rm -f "$TEMP_DIR"/chunk_* $(cat "$TEMP_DIR/chunk_list.txt")
fi

NORMALIZED_COUNT=$(wc -l < "$TEMP_DIR/normalized.txt")
echo -e "${GREEN}  ✓ Normalized edges: $NORMALIZED_COUNT${NC}"
echo "  ✓ Self-loops removed: $((TOTAL_LINES - NORMALIZED_COUNT - 1))"

# ========================================
# STEP 2: Sort and deduplicate
# ========================================
echo -e "${YELLOW}Step 2: Sorting and removing duplicates...${NC}"

if [ "$USE_MEMORY_SORT" = true ]; then
    echo "  Using memory-based sorting..."
    
    # Sort by first two columns, then by weight (descending)
    # Keep only the highest weight for each edge pair
    sort -k1,1 -k2,2 -k3,3nr "$TEMP_DIR/normalized.txt" | \
        awk '{
            edge = $1 "_" $2
            if (!(edge in seen)) {
                print $1, $2, $3
                seen[edge] = 1
            }
        }' > "$OUTPUT_FILE"
    
else
    echo "  Using disk-based sorting with parallel merge..."
    
    # Use GNU sort with parallel option
    export LC_ALL=C  # Faster sorting
    
    sort --parallel="$NUM_JOBS" \
         --buffer-size=2G \
         --temporary-directory="$TEMP_DIR" \
         -k1,1 -k2,2 -k3,3nr \
         "$TEMP_DIR/normalized.txt" | \
        awk '{
            edge = $1 "_" $2
            if (!(edge in seen)) {
                print $1, $2, $3
                seen[edge] = 1
            }
        }' > "$OUTPUT_FILE"
fi

# ========================================
# STEP 3: Statistics and validation
# ========================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

FINAL_COUNT=$(wc -l < "$OUTPUT_FILE")
REMOVED_DUPLICATES=$((NORMALIZED_COUNT - FINAL_COUNT))

echo ""
echo -e "${GREEN}===================================="
echo " EDGELIST GENERATION COMPLETE"
echo "====================================${NC}"
echo " Processing time: ${ELAPSED}s"
echo " Input edges:     $((TOTAL_LINES - 1))"
echo " Self-loops removed: $((TOTAL_LINES - NORMALIZED_COUNT - 1))"
echo " Duplicates removed: $REMOVED_DUPLICATES"
echo " Final edges:     $FINAL_COUNT"
echo " Reduction:       $(echo "scale=1; (1 - $FINAL_COUNT / ($TOTAL_LINES - 1)) * 100" | bc)%"
echo -e "${GREEN}====================================${NC}"
echo " Output saved to:"
echo " $OUTPUT_FILE"
echo -e "${GREEN}====================================${NC}"

# Additional analysis
if [ "$FINAL_COUNT" -gt 0 ]; then
    echo ""
    echo "Edge weight statistics:"
    awk '{sum+=$3; if(NR==1||$3>max)max=$3; if(NR==1||$3<min)min=$3} 
         END {printf "  Min weight:  %.1f\n  Max weight:  %.1f\n  Mean weight: %.1f\n", min, max, sum/NR}' \
         "$OUTPUT_FILE"
    
    # Check connectivity
    echo ""
    echo "Network connectivity:"
    UNIQUE_NODES=$(awk '{print $1; print $2}' "$OUTPUT_FILE" | sort -u | wc -l)
    AVG_DEGREE=$(echo "scale=1; $FINAL_COUNT * 2 / $UNIQUE_NODES" | bc)
    
    echo "  Unique nodes: $UNIQUE_NODES"
    echo "  Average degree: $AVG_DEGREE"
    
    # Find highest degree node
    echo -n "  Max degree node: "
    awk '{print $1; print $2}' "$OUTPUT_FILE" | \
        sort | uniq -c | sort -rn | head -1 | \
        awk '{printf "%s (degree: %d)\n", $2, $1}'
    
    # Warning for problematic clustering
    if (( $(echo "$AVG_DEGREE > 100" | bc -l) )); then
        echo ""
        echo -e "${YELLOW}⚠ WARNING: High average degree detected!${NC}"
        echo "  This may cause large clusters in MCL."
        echo "  Consider:"
        echo "    1. Stricter BLAST filtering"
        echo "    2. Higher MCL inflation (2.5-3.0)"
        echo "    3. Edge weight threshold filtering"
    fi
fi

exit 0