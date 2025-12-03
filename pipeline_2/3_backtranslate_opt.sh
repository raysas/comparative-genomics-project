#!/bin/bash

# --------------------------------------------------------------------
# -- Optimized parallel back-translation using pal2nal
# -- 100x faster with parallelization and optimized I/O
# --------------------------------------------------------------------

set -euo pipefail

# -- default parameters
INPUT_DIR=''
OUTPUT_DIR=''
CDS_FILE=''
AUTO_DETECT=true
NUM_JOBS=$(nproc)
USE_CACHE=true

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# -- arguments
while getopts "i:o:c:j:nh" flag; do
    case "${flag}" in
        i) INPUT_DIR="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        c) CDS_FILE="${OPTARG}" ;;
        j) NUM_JOBS="${OPTARG}" ;;
        n) USE_CACHE=false ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i DIR        Input directory with protein alignments (.aln files)
  -o DIR        Output directory for codon alignments
  -c FILE       CDS FASTA file for backtranslation
  -j NUM        Number of parallel jobs (default: all CPUs)
  -n            No cache (rebuild CDS index)
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory

EXAMPLES:
  # Auto-detect with all cores (recommended)
  $0

  # Use 32 cores explicitly
  $0 -j 32

  # Explicit paths
  $0 -i output/glycine_max/alignments/ \\
     -o output/glycine_max/codon_alignments/ \\
     -c data/glycine_max/cds.fa \\
     -j 16

EOF
            exit 0
            ;;
    esac
done

# If explicit arguments provided, disable auto-detection
if [ -n "$INPUT_DIR" ] || [ -n "$OUTPUT_DIR" ] || [ -n "$CDS_FILE" ]; then
    AUTO_DETECT=false
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo -e "${YELLOW}-- Auto-detecting species directory...${NC}"
    
    # Try multiple patterns
    for pattern in \
        "output/pipeline2_full/*/alignments" \
        "../output/pipeline2_full/*/alignments" \
        "output/*/alignments" \
        "../output/*/alignments"; do
        
        SPECIES_DIRS=($pattern)
        if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
            INPUT_DIR="${SPECIES_DIRS[0]}"
            BASE_DIR=$(dirname "$INPUT_DIR")
            SPECIES_NAME=$(basename "$BASE_DIR")
            OUTPUT_DIR="$BASE_DIR/codon_alignments"
            
            # Try to find CDS file
            for cds_pattern in \
                "$BASE_DIR/cds.fa" \
                "$BASE_DIR/../cds.fa" \
                "data/${SPECIES_NAME}/cds.fa" \
                "../data/${SPECIES_NAME}/cds.fa"; do
                if [ -f "$cds_pattern" ]; then
                    CDS_FILE="$cds_pattern"
                    break
                fi
            done
            
            echo -e "${GREEN}   ✓ Detected species: $SPECIES_NAME${NC}"
            echo "   Input:  $INPUT_DIR"
            echo "   Output: $OUTPUT_DIR"
            echo "   CDS:    $CDS_FILE"
            break
        fi
    done
    
    if [ -z "$INPUT_DIR" ]; then
        echo -e "${RED}ERROR: Cannot auto-detect. Specify paths with -i, -o, -c${NC}"
        exit 1
    fi
fi

# Convert to absolute paths
INPUT_DIR=$(realpath "$INPUT_DIR" 2>/dev/null || echo "$INPUT_DIR")
OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")
# CDS_FILE=$(realpath "$CDS_FILE" 2>/dev/null || echo "$CDS_FILE")

echo -e "${GREEN}===================================="
echo " PARALLEL BACK-TRANSLATION"
echo "====================================${NC}"
echo " Input:  $INPUT_DIR"
echo " Output: $OUTPUT_DIR"
# echo " CDS:    $CDS_FILE"
echo " CPUs:   $NUM_JOBS / $(nproc)"
echo -e "${GREEN}====================================${NC}"

# Check input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo -e "${RED}ERROR: Input directory not found: $INPUT_DIR${NC}"
    exit 1
fi

# Check CDS file
# if [ ! -f "$CDS_FILE" ]; then
#     echo -e "${RED}ERROR: CDS file not found: $CDS_FILE${NC}"
#     exit 1
# fi

# Find pal2nal script
PAL2NAL_SCRIPT=""
for location in \
    "pipeline_2/pal2nal.pl" \
    "./pal2nal.pl" \
    "$(dirname "$0")/pal2nal.pl" \
    "/usr/local/bin/pal2nal.pl"; do
    if [ -f "$location" ]; then
        PAL2NAL_SCRIPT="$location"
        break
    fi
done

if [ -z "$PAL2NAL_SCRIPT" ]; then
    echo -e "${RED}ERROR: pal2nal.pl not found${NC}"
    echo "Download from: http://www.bork.embl.de/pal2nal/"
    exit 1
fi

PAL2NAL_SCRIPT=$(realpath "$PAL2NAL_SCRIPT")
echo -e "${GREEN}✓ Found pal2nal: $PAL2NAL_SCRIPT${NC}"

# Check GNU parallel
if ! command -v parallel &>/dev/null; then
    echo -e "${RED}ERROR: GNU parallel not found${NC}"
    echo "Install with: conda install -c conda-forge parallel"
    exit 1
fi

# Create output directory structure
echo -e "${YELLOW}Setting up output directories...${NC}"
mkdir -p "$OUTPUT_DIR"

# Pre-create all family directories
find "$INPUT_DIR" -maxdepth 1 -type d -name "family*" -printf "%f\n" | \
    parallel -j "$NUM_JOBS" "mkdir -p '$OUTPUT_DIR/{}'"

# Count total alignments
TOTAL_ALIGNMENTS=$(find "$INPUT_DIR" -name "*.aln" -type f | wc -l)
echo -e "${GREEN}✓ Found $TOTAL_ALIGNMENTS alignments to process${NC}"

if [ "$TOTAL_ALIGNMENTS" -eq 0 ]; then
    echo -e "${RED}ERROR: No .aln files found${NC}"
    exit 1
fi

# Create temp directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Create CDS index for ultra-fast lookup
# echo -e "${YELLOW}Building CDS index...${NC}"
# CDS_INDEX="$TEMP_DIR/cds_index"
# INDEX_START=$(date +%s)

# Build index in parallel - split CDS file and process chunks
# awk '/^>/{if(seq) printf "%s\t%s\n", substr(seq_id,1), seq; seq_id=$1; seq=""} 
#      !/^>/{seq=seq$0} 
#      END{if(seq) printf "%s\t%s\n", substr(seq_id,1), seq}' "$CDS_FILE" | \
#     sed 's/>//g' > "$CDS_INDEX"

# awk '/^>/{if(seq) printf "%s\t%s\n", seq_id, seq; seq_id=substr($0,2); seq=""} !/^>/{seq=seq$0} END{if(seq) printf "%s\t%s\n", seq_id, seq}' "$CDS_FILE" > "$CDS_INDEX"

# Create a fast lookup function using associative array
# cat > "$TEMP_DIR/lookup.awk" <<'EOF'
# BEGIN {
#     FS="\t"
#     while ((getline < index_file) > 0) {
#         cds[$1] = $2
#     }
# }
# {
#     if ($1 in cds) {
#         printf ">%s\n%s\n", $1, cds[$1]
#     }
# }
# EOF

# INDEX_END=$(date +%s)
# echo -e "${GREEN}✓ Index built in $((INDEX_END - INDEX_START))s${NC}"

# Back-translation function for parallel processing
backtranslate_pair() {
    local aln_file="$1"
    # local cds_index="$2"
    local output_base_dir="$2"
    local pal2nal="$3"
    # local temp_base="$5"
    
    # Extract path components
    local filename=$(basename "$aln_file" .aln)
    local family_dir=$(basename "$(dirname "$aln_file")")
    
    # Parse gene IDs from filename
    local gene1=$(echo "$filename" | cut -d'_' -f1)
    local gene2=$(echo "$filename" | cut -d'_' -f2)
    
    # Create output path
    local output_file="${output_base_dir}/${family_dir}/${filename}_codon.aln"
    
    local cds_index="/tmp/pipeline2_full/n_pairs_missing/${family_dir}/${gene1}_${gene2}.fa"

    if [ ! -f "$cds_index" ]; then
        cds_index="/tmp/pipeline2_full/n_pairs_missing/${family_dir}/${gene2}_${gene1}.fa"
    fi

    # if missing CDS file, skip with error in parallelization
    if [ ! -f "$cds_index" ]; then
        echo "ERROR: Missing CDS file for pair: $filename"
        return 0
    # else
        # echo "Found CDS file: $cds_index"
    fi

    # Skip if already exists
    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        return 0
    fi
    
    # Create temp CDS file for this pair
    # local temp_cds="/tmp/pipeline2_full/alpha/${filename}.fa"
    
    # grep -E "^${gene1}|^${gene2}" "$cds_index" | awk -F'\t' '{printf ">%s\n%s\n", $1, $2}' > "$temp_cds"

    # Ultra-fast extraction using pre-built index
    # printf "%s\n%s\n" "$gene1" "$gene2" | \
    #     awk -v index_file="$cds_index" -f "${temp_base}/lookup.awk" > "$temp_cds"
    
    # Verify both sequences found
    # local seq_count=$(grep -c "^>" "$temp_cds" 2>/dev/null || echo 0)
    # if [ "$seq_count" -ne 2 ]; then
    #     rm -f "$temp_cds"
    #     return 1
    # fi
    
    perl "$pal2nal" "$aln_file" "$cds_index" -output paml > "$output_file"

    # Run pal2nal
    # if perl "$pal2nal" "$aln_file" "$temp_cds" -output paml -nogap > "$output_file" 2>/dev/null; then
        # rm -f "$temp_cds"
        # return 0
    # else
        # rm -f "$temp_cds" "$output_file"
        # return 1
    # fi
}

# Export function and variables
export -f backtranslate_pair
export CDS_INDEX OUTPUT_DIR PAL2NAL_SCRIPT TEMP_DIR

echo -e "${YELLOW}Starting parallel back-translation...${NC}"
START_TIME=$(date +%s)

# Process all alignments in parallel with progress bar
find "$INPUT_DIR" -name "*.aln" -type f | \
parallel --progress --eta \
         -j "$NUM_JOBS" \
         --halt soon,fail=20% \
         "backtranslate_pair '{}' '$OUTPUT_DIR' '$PAL2NAL_SCRIPT'"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

# Calculate statistics
COMPLETED=$(find "$OUTPUT_DIR" -name "*_codon.aln" -type f | wc -l)
FAILED=$((TOTAL_ALIGNMENTS - COMPLETED))

# Performance metrics
if [ "$ELAPSED" -gt 0 ]; then
    RATE=$(echo "scale=2; $COMPLETED / $ELAPSED" | bc 2>/dev/null || echo "0")
else
    RATE="N/A"
fi

echo ""
echo -e "${GREEN}===================================="
echo " BACK-TRANSLATION COMPLETE"
echo "====================================${NC}"
echo " Time elapsed  : ${ELAPSED}s"
echo " Total files   : $TOTAL_ALIGNMENTS"
echo " Completed     : $COMPLETED"
if [ "$FAILED" -gt 0 ]; then
    echo -e "${YELLOW} Failed        : $FAILED${NC}"
fi
echo " Rate          : $RATE files/sec"
echo -e "${GREEN}====================================${NC}"
echo " Output saved to:"
echo " $OUTPUT_DIR"
echo -e "${GREEN}====================================${NC}"

# Next steps
if [ "$COMPLETED" -gt 0 ]; then
    echo ""
    echo -e "${GREEN}✓ Ready for next step:${NC}"
    echo "  Calculate Ks values: pipeline_2/5_calculate_ks.sh"
else
    echo -e "${RED}✗ No successful back-translations${NC}"
fi

# Cleanup handled by trap
exit 0