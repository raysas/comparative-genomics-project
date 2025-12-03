#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    Aligns pairwise protein sequences using clustalw2 
#    (input from 1_prepare_pairs.sh output)
# -- Usage:
#    bash ./pipeline_2/3_align_proteins.sh [-i INPUT_DIR] [-o OUTPUT_DIR] [-h]
# -- default (without params) equivalent to:
#    Auto-detects from species directories
# --------------------------------------------------------------------

# -- default parameters
INPUT_DIR=''
OUTPUT_DIR=''
AUTO_DETECT=true

# -- arguments
while getopts "i:o:h" flag; do
    case "${flag}" in
        i) INPUT_DIR="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i DIR        Input directory with family pairs (families/)
  -o DIR        Output directory for alignments
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
    Input:  output/{species}/families/
    Output: output/{species}/alignments/

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i output/glycine_max/families/ \\
     -o output/glycine_max/alignments/

EOF
            exit 0
            ;;
        *)
            echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
    esac
done

# If explicit arguments provided, disable auto-detection
if [ -n "$INPUT_DIR" ] || [ -n "$OUTPUT_DIR" ]; then
    AUTO_DETECT=false
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    SPECIES_DIRS=(output/*/families)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        INPUT_DIR="${SPECIES_DIRS[0]}"
        SPECIES_NAME=$(basename $(dirname "$INPUT_DIR"))
        OUTPUT_DIR="output/${SPECIES_NAME}/alignments"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Input:  $INPUT_DIR"
        echo "   Output: $OUTPUT_DIR"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species family directories found"
        echo "       Specify paths explicitly with -i and -o"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No family pairs found in output/*/"
        echo "       Run pipeline_2/1_prepare_pairs.sh first"
        exit 1
    fi
fi

# Setup logging
LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo "-- Parameters:"
echo "   INPUT DIR  : $INPUT_DIR"
echo "   OUTPUT DIR : $OUTPUT_DIR"

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    echo "       Run pipeline_2/1_prepare_pairs.sh first"
    exit 1
fi

# Check if clustalw is available (prefer clustalw2, fallback to clustalw)
CLUSTALW_CMD=""
if command -v clustalw2 &> /dev/null; then
    CLUSTALW_CMD="clustalw2"
elif command -v clustalw &> /dev/null; then
    CLUSTALW_CMD="clustalw"
else
    echo "ERROR: ClustalW not found. Please install:"
    echo "       conda install -c bioconda clustalw"
    echo "   or  apt-get install clustalw"
    exit 1
fi

echo "   Using alignment tool: $CLUSTALW_CMD"

echo "-- Aligning protein pairs with ClustalW2..."

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count families and pairs for progress
FAMILY_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -name "family*" -type d | wc -l)
TOTAL_PAIRS=$(find "$INPUT_DIR" -name "*.fa" | wc -l)

echo "   Found $FAMILY_COUNT families with $TOTAL_PAIRS protein pairs"

if [ $TOTAL_PAIRS -eq 0 ]; then
    echo "ERROR: No protein pairs found in $INPUT_DIR"
    echo "       Run pipeline_2/1_prepare_pairs.sh first"
    exit 1
fi

START_TIME=$(date +%s)
PROCESSED_PAIRS=0

# Process all family directories
for fam_dir in "$INPUT_DIR"/family*; do
    if [ ! -d "$fam_dir" ]; then
        continue
    fi
    
    fam_name=$(basename "$fam_dir")
    fam_out="$OUTPUT_DIR/$fam_name"
    mkdir -p "$fam_out"

    echo "▶ Processing $fam_name"

    # Count pairs in this family
    FAMILY_PAIRS=$(find "$fam_dir" -name "*.fa" | wc -l)
    
    if [ $FAMILY_PAIRS -eq 0 ]; then
        echo "   No pairs found in $fam_name"
        continue
    fi

    # Process all pairwise FASTA files in this family
    for pair in "$fam_dir"/*.fa; do
        if [ ! -f "$pair" ]; then
            continue
        fi
        
        pair_base=$(basename "$pair" .fa)
        out_file="$fam_out/${pair_base}.aln"
        
        # Skip if alignment already exists
        if [ -f "$out_file" ]; then
            echo "   ✓ $pair_base (already aligned)"
            ((PROCESSED_PAIRS++))
            continue
        fi
        
        echo "   ▶ Aligning $pair_base"
        
        # Run ClustalW alignment
        if $CLUSTALW_CMD -quiet -align -infile="$pair" -outfile="$out_file" >/dev/null 2>&1; then
            echo "   ✓ $pair_base -> ${pair_base}.aln"
        else
            echo "   ✗ Failed to align $pair_base"
        fi
        
        ((PROCESSED_PAIRS++))
        
        # Progress update every 100 pairs
        if [ $((PROCESSED_PAIRS % 100)) -eq 0 ]; then
            ELAPSED=$(($(date +%s) - START_TIME))
            echo "   Progress: $PROCESSED_PAIRS/$TOTAL_PAIRS pairs (${ELAPSED}s elapsed)"
        fi
    done
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "-- Protein alignment completed successfully"
echo "   Total pairs processed: $PROCESSED_PAIRS"
echo "   Processing time: ${ELAPSED}s"
echo "   Output directory: $OUTPUT_DIR"

echo ""
echo "-- Next steps:"
echo "   1. Back-translate alignments (pipeline_2/4_backtranslate.sh)"
echo "   2. Calculate Ks values (pipeline_2/6_run_yn00.sh)"

exit 0