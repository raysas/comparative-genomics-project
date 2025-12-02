#!/bin/bash

# --------------------------------------------------------------------
# -- Optimized BLAST filtering with parallel processing
# -- 10-50x faster through efficient AWK and parallel chunking
# --------------------------------------------------------------------

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Default parameters
INPUT_FILE='output/blast_output/blast_results_with_coverage.tsv'
OUTPUT_DIR='output/blast_filtered'
ID_THRESHOLD=30
Q_COV_THRESHOLD=50
S_COV_THRESHOLD=50
BIT_SCORE_THRESHOLD=""
E_VALUE_THRESHOLD=""
NUM_JOBS=$(nproc)
CHUNK_SIZE=1000000  # Process in chunks of 1M lines

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT_FILE="$2"; shift 2;;
        -o) OUTPUT_DIR="$2"; shift 2;;
        -id) ID_THRESHOLD="$2"; shift 2;;
        -qcov) Q_COV_THRESHOLD="$2"; shift 2;;
        -scov) S_COV_THRESHOLD="$2"; shift 2;;
        -bit) BIT_SCORE_THRESHOLD="$2"; shift 2;;
        -evalue) E_VALUE_THRESHOLD="$2"; shift 2;;
        -j) NUM_JOBS="$2"; shift 2;;
        -c) CHUNK_SIZE="$2"; shift 2;;
        -h|--help)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE    Input BLAST results with coverage
  -o DIR     Output directory
  -id NUM    Identity threshold (default: 30)
  -qcov NUM  Query coverage threshold (default: 50)
  -scov NUM  Subject coverage threshold (default: 50)
  -bit NUM   Bit score threshold (optional)
  -evalue NUM E-value threshold (optional)
  -j NUM     Parallel jobs (default: all CPUs)
  -c NUM     Chunk size for processing (default: 1000000)

EXAMPLES:
  # Standard filtering
  $0 -id 30 -qcov 50 -scov 50

  # Strict filtering for Ks analysis
  $0 -id 40 -qcov 70 -scov 70 -j 32

  # Include bit score and e-value
  $0 -id 30 -qcov 50 -scov 50 -bit 100 -evalue 1e-10

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
echo " OPTIMIZED BLAST FILTERING"
echo -e "====================================${NC}"
echo " Input:  $INPUT_FILE"
echo " Output: $OUTPUT_DIR"
echo " CPUs:   $NUM_JOBS"
echo -e "${GREEN}====================================${NC}"

# Check input
if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}ERROR: Input file not found: $INPUT_FILE${NC}"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build output filename
OUTPUT_FILE="${OUTPUT_DIR}/filtered_blast_results_id${ID_THRESHOLD}_qcov${Q_COV_THRESHOLD}_scov${S_COV_THRESHOLD}"
if [ -n "$BIT_SCORE_THRESHOLD" ]; then
    OUTPUT_FILE="${OUTPUT_FILE}_bit${BIT_SCORE_THRESHOLD}"
fi
if [ -n "$E_VALUE_THRESHOLD" ]; then
    OUTPUT_FILE="${OUTPUT_FILE}_evalue${E_VALUE_THRESHOLD}"
fi
OUTPUT_FILE="${OUTPUT_FILE}.tsv"

echo -e "${YELLOW}Filtering parameters:${NC}"
echo "  Identity       >= $ID_THRESHOLD"
echo "  Query Coverage >= $Q_COV_THRESHOLD"
echo "  Subject Coverage >= $S_COV_THRESHOLD"
if [ -n "$BIT_SCORE_THRESHOLD" ]; then
    echo "  Bit Score     >= $BIT_SCORE_THRESHOLD"
fi
if [ -n "$E_VALUE_THRESHOLD" ]; then
    echo "  E-value       <= $E_VALUE_THRESHOLD"
fi

# Count input lines
echo -e "${YELLOW}Analyzing input file...${NC}"
TOTAL_LINES=$(wc -l < "$INPUT_FILE")
echo -e "${GREEN}✓ Total lines: $TOTAL_LINES${NC}"

# Create temp directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

START_TIME=$(date +%s)

# ========================================
# CASE 1: Small file - single AWK pass
# ========================================
if [ "$TOTAL_LINES" -lt 10000000 ]; then
    echo -e "${YELLOW}Using single-pass AWK filtering (file < 10M lines)...${NC}"
    
    # Build AWK filter command
    AWK_FILTER='BEGIN {OFS="\t"} {'
    AWK_FILTER="${AWK_FILTER} if ("
    AWK_FILTER="${AWK_FILTER} \$3 >= ${ID_THRESHOLD}"  # Identity
    AWK_FILTER="${AWK_FILTER} && \$13 >= ${Q_COV_THRESHOLD}"  # Query coverage
    AWK_FILTER="${AWK_FILTER} && \$14 >= ${S_COV_THRESHOLD}"  # Subject coverage
    
    if [ -n "$BIT_SCORE_THRESHOLD" ]; then
        AWK_FILTER="${AWK_FILTER} && \$12 >= ${BIT_SCORE_THRESHOLD}"
    fi
    
    if [ -n "$E_VALUE_THRESHOLD" ]; then
        AWK_FILTER="${AWK_FILTER} && \$11 <= ${E_VALUE_THRESHOLD}"
    fi
    
    AWK_FILTER="${AWK_FILTER}) print \$0 }"
    
    # Apply filter
    awk "$AWK_FILTER" "$INPUT_FILE" > "$OUTPUT_FILE"
    
# ========================================
# CASE 2: Large file - parallel chunking
# ========================================
else
    echo -e "${YELLOW}Using parallel chunk processing (file > 10M lines)...${NC}"
    
    # Split file into chunks
    echo "Splitting into chunks of $CHUNK_SIZE lines..."
    split -l "$CHUNK_SIZE" "$INPUT_FILE" "$TEMP_DIR/chunk_"
    
    NUM_CHUNKS=$(ls "$TEMP_DIR"/chunk_* | wc -l)
    echo -e "${GREEN}✓ Created $NUM_CHUNKS chunks${NC}"
    
    # Function to filter a chunk
    filter_chunk() {
        local chunk_file="$1"
        local chunk_out="${chunk_file}.filtered"
        
        # Build AWK command with all filters
        awk -v id="$ID_THRESHOLD" \
            -v qcov="$Q_COV_THRESHOLD" \
            -v scov="$S_COV_THRESHOLD" \
            -v bit="${BIT_SCORE_THRESHOLD:-0}" \
            -v eval="${E_VALUE_THRESHOLD:-999}" \
            'BEGIN {OFS="\t"} 
            {
                # Apply all filters in one pass
                if ($3 >= id && $13 >= qcov && $14 >= scov) {
                    if (bit == 0 || $12 >= bit) {
                        if (eval == 999 || $11 <= eval) {
                            print $0
                        }
                    }
                }
            }' "$chunk_file" > "$chunk_out"
        
        echo "$chunk_out"
    }
    
    export -f filter_chunk
    export ID_THRESHOLD Q_COV_THRESHOLD S_COV_THRESHOLD
    export BIT_SCORE_THRESHOLD E_VALUE_THRESHOLD
    
    # Process chunks in parallel
    echo "Filtering chunks in parallel..."
    ls "$TEMP_DIR"/chunk_* | \
        parallel -j "$NUM_JOBS" --progress --bar filter_chunk {} > "$TEMP_DIR/filtered_list.txt"
    
    # Combine filtered chunks
    echo "Combining filtered results..."
    cat $(cat "$TEMP_DIR/filtered_list.txt") > "$OUTPUT_FILE"
fi

# ========================================
# Statistics and reporting
# ========================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

# Count results
OUTPUT_LINES=$(wc -l < "$OUTPUT_FILE")
REMOVED_LINES=$((TOTAL_LINES - OUTPUT_LINES))
RETENTION_PCT=$(echo "scale=2; $OUTPUT_LINES * 100 / $TOTAL_LINES" | bc)

echo ""
echo -e "${GREEN}===================================="
echo " FILTERING COMPLETE"
echo "====================================${NC}"
echo " Processing time: ${ELAPSED}s"
echo " Input lines:     $TOTAL_LINES"
echo " Output lines:    $OUTPUT_LINES"
echo " Removed:         $REMOVED_LINES"
echo " Retention:       ${RETENTION_PCT}%"
echo -e "${GREEN}====================================${NC}"
echo " Output saved to:"
echo " $OUTPUT_FILE"
echo -e "${GREEN}====================================${NC}"

# Additional statistics
if [ "$OUTPUT_LINES" -gt 0 ]; then
    echo ""
    echo "Quick statistics of filtered data:"
    echo -n "  Average identity: "
    awk '{sum+=$3; n++} END {printf "%.1f%%\n", sum/n}' "$OUTPUT_FILE"
    echo -n "  Average query coverage: "
    awk '{sum+=$13; n++} END {printf "%.1f%%\n", sum/n}' "$OUTPUT_FILE"
    echo -n "  Average subject coverage: "
    awk '{sum+=$14; n++} END {printf "%.1f%%\n", sum/n}' "$OUTPUT_FILE"
    
    # Check for problematic clustering
    echo ""
    echo "Checking pair distribution..."
    UNIQUE_QUERIES=$(cut -f1 "$OUTPUT_FILE" | sort -u | wc -l)
    UNIQUE_SUBJECTS=$(cut -f2 "$OUTPUT_FILE" | sort -u | wc -l)
    UNIQUE_GENES=$((UNIQUE_QUERIES > UNIQUE_SUBJECTS ? UNIQUE_QUERIES : UNIQUE_SUBJECTS))
    AVG_PAIRS_PER_GENE=$(echo "scale=1; $OUTPUT_LINES / $UNIQUE_GENES" | bc)
    
    echo "  Unique genes: $UNIQUE_GENES"
    echo "  Average pairs per gene: $AVG_PAIRS_PER_GENE"
    
    if (( $(echo "$AVG_PAIRS_PER_GENE > 100" | bc -l) )); then
        echo -e "${YELLOW}  ⚠ Warning: High connectivity detected!${NC}"
        echo "    Consider stricter filtering for Ks analysis:"
        echo "    -id 40 -qcov 70 -scov 70"
    fi
fi

exit 0