#!/bin/bash

# --------------------------------------------------------------------
# -- Robust parallel MAFFT alignment with proper path handling
# --------------------------------------------------------------------

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# -- default parameters
INPUT_DIR=''
OUTPUT_DIR=''
NUM_JOBS=$(nproc)
DEBUG=false

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# -- arguments
while getopts "i:o:j:dh" flag; do
    case "${flag}" in
        i) INPUT_DIR="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        j) NUM_JOBS="${OPTARG}" ;;
        d) DEBUG=true ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i DIR   Input directory with family*/pair.fa structure
  -o DIR   Output directory for alignments  
  -j NUM   Number of parallel jobs (default: $(nproc))
  -d       Debug mode (verbose output)
  -h       Show this help

EXAMPLES:
  # Auto-detect paths
  $0

  # Explicit paths
  $0 -i ../output/pipeline2_full/glycine_max/pairs \\
     -o ../output/pipeline2_full/glycine_max/alignments

  # Debug mode with 8 cores
  $0 -d -j 8

EOF
            exit 0
            ;;
        *)
            echo "Invalid option. Use -h for help."
            exit 1
            ;;
    esac
done

# Function to print colored messages
print_msg() {
    local color=$1
    shift
    echo -e "${color}$*${NC}"
}

# Auto-detect paths if not provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    print_msg "$YELLOW" "Auto-detecting directories..."
    
    # Try common patterns
    for pattern in \
        "output/pipeline2_full/*/pairs" \
        "output/pipeline2_full/*/families" \
        "../output/pipeline2_full/*/pairs" \
        "../output/pipeline2_full/*/families" \
        "output/*/pairs" \
        "output/*/families"; do
        
        for dir in $pattern; do
            if [ -d "$dir" ]; then
                INPUT_DIR="$dir"
                # Replace pairs/families with alignments
                OUTPUT_DIR="${dir%/*}/alignments"
                print_msg "$GREEN" "✓ Found: $INPUT_DIR"
                break 2
            fi
        done
    done
    
    if [ -z "$INPUT_DIR" ]; then
        print_msg "$RED" "ERROR: Cannot find input directory"
        echo "Please specify with -i flag"
        exit 1
    fi
fi

# Resolve to absolute paths
INPUT_DIR=$(cd "$INPUT_DIR" 2>/dev/null && pwd || realpath "$INPUT_DIR" 2>/dev/null || echo "$INPUT_DIR")
OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")

print_msg "$GREEN" "===================================="
print_msg "$GREEN" " PARALLEL MAFFT ALIGNMENT"
print_msg "$GREEN" "===================================="
echo " Input:  $INPUT_DIR"
echo " Output: $OUTPUT_DIR"
echo " CPUs:   $NUM_JOBS / $(nproc) available"
print_msg "$GREEN" "===================================="

# Verify dependencies
check_tool() {
    if ! command -v "$1" &>/dev/null; then
        print_msg "$RED" "ERROR: $1 not found"
        echo "Install with: $2"
        exit 1
    fi
}

check_tool "mafft" "conda install -c bioconda mafft"
check_tool "parallel" "conda install -c conda-forge parallel"

# Verify input directory
if [ ! -d "$INPUT_DIR" ]; then
    print_msg "$RED" "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Count input files and create directory structure
print_msg "$YELLOW" "Analyzing input structure..."

# Get list of family directories
FAMILY_DIRS=$(find "$INPUT_DIR" -maxdepth 1 -type d -name "family*" | sort)
FAMILY_COUNT=$(echo "$FAMILY_DIRS" | grep -c . || echo 0)

if [ "$FAMILY_COUNT" -eq 0 ]; then
    print_msg "$RED" "ERROR: No family* directories found in $INPUT_DIR"
    exit 1
fi

# Create output directory structure
mkdir -p "$OUTPUT_DIR"
echo "$FAMILY_DIRS" | while read fam_dir; do
    [ -z "$fam_dir" ] && continue
    family_name=$(basename "$fam_dir")
    mkdir -p "$OUTPUT_DIR/$family_name"
done

# Count total files
TOTAL_FILES=$(find "$INPUT_DIR" -name "*.fa" -type f | wc -l)
print_msg "$GREEN" "✓ Found $FAMILY_COUNT families with $TOTAL_FILES pairs"

if [ "$TOTAL_FILES" -eq 0 ]; then
    print_msg "$RED" "ERROR: No .fa files found"
    exit 1
fi

# Create alignment function
align_file() {
    local input_file="$1"
    
    # Debug output
    [ "$DEBUG" = true ] && echo "Processing: $input_file"
    
    # Extract components from path
    local filename=$(basename "$input_file")
    local family_dir=$(basename "$(dirname "$input_file")")
    local base_name="${filename%.fa}"
    
    # Build output path
    local output_file="$OUTPUT_DIR/${family_dir}/${base_name}.aln"
    
    # Skip if already exists and is non-empty
    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        [ "$DEBUG" = true ] && echo "  Skip (exists): $output_file"
        return 0
    fi
    
    # Create temporary file in the same directory as output
    local temp_file="${output_file}.tmp.$$.$RANDOM"
    
    # Run MAFFT with speed optimizations
    if mafft --amino --retree 1 --maxiterate 0 --thread 1 --quiet \
             "$input_file" > "$temp_file" 2>/dev/null; then
        mv "$temp_file" "$output_file"
        [ "$DEBUG" = true ] && echo "  Created: $output_file"
        return 0
    else
        rm -f "$temp_file"
        [ "$DEBUG" = true ] && echo "  Failed: $input_file"
        return 1
    fi
}

# Export for parallel
export -f align_file
export OUTPUT_DIR DEBUG

# Run alignment
print_msg "$YELLOW" "Starting alignment process..."
START_TIME=$(date +%s)

# Create a progress bar and run in parallel
if [ "$DEBUG" = true ]; then
    # Debug mode - verbose output
    find "$INPUT_DIR" -name "*.fa" -type f | \
        parallel -j "$NUM_JOBS" \
                 --halt soon,fail=20% \
                 align_file
else
    # Normal mode - with progress bar
    find "$INPUT_DIR" -name "*.fa" -type f | \
        parallel -j "$NUM_JOBS" \
                 --progress --colsep '\t' \
                 --halt soon,fail=20% \
                 align_file
fi

# Calculate statistics
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

# Count completed alignments
COMPLETED=$(find "$OUTPUT_DIR" -name "*.aln" -type f | wc -l)
FAILED=$((TOTAL_FILES - COMPLETED))

# Calculate rate
if [ "$ELAPSED" -gt 0 ]; then
    RATE=$(echo "scale=2; $COMPLETED / $ELAPSED" | bc 2>/dev/null || echo "0")
else
    RATE="N/A"
fi

# Final report
echo ""
print_msg "$GREEN" "===================================="
print_msg "$GREEN" " ALIGNMENT COMPLETE"
print_msg "$GREEN" "===================================="
echo " Time elapsed : ${ELAPSED}s"
echo " Total files  : $TOTAL_FILES"
echo " Completed    : $COMPLETED"
if [ "$FAILED" -gt 0 ]; then
    print_msg "$YELLOW" " Failed       : $FAILED"
fi
echo " Rate         : $RATE alignments/sec"
print_msg "$GREEN" "===================================="
echo " Output saved to:"
echo " $OUTPUT_DIR"
print_msg "$GREEN" "===================================="

# Provide next steps
if [ "$FAILED" -eq 0 ]; then
    print_msg "$GREEN" "✓ All alignments completed successfully!"
else
    print_msg "$YELLOW" "⚠ Some alignments failed. Re-run to retry."
fi

exit 0