#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    Consolidates PAML yn00 results into a summary table with QC filters
#    (parses individual yn00 outputs and creates final Ks dataset)
# -- Usage:
#    bash ./pipeline_2/6_consolidate_results.sh [-i INPUT_DIR] [-o OUTPUT_FILE] [-h]
# -- default (without params) equivalent to:
#    Auto-detects from species directories
# --------------------------------------------------------------------

# -- default parameters
INPUT_DIR=''
OUTPUT_FILE=''
AUTO_DETECT=true

# Quality control parameters
MAX_KS=5.0
MAX_KA=5.0
MAX_KAKS=5.0
MIN_ALIGNMENT_LENGTH=50

# -- arguments
while getopts "i:o:h" flag; do
    case "${flag}" in
        i) INPUT_DIR="${OPTARG}" ;;
        o) OUTPUT_FILE="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i DIR        Input directory with yn00 results
  -o FILE       Output summary table file (.tsv)
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
    Input:  output/{species}/ks_results/
    Output: output/{species}/ks_summary.tsv

QUALITY FILTERS:
  - Maximum Ks: $MAX_KS (filters saturated pairs)
  - Maximum Ka: $MAX_KA (filters unrealistic values)
  - Maximum Ka/Ks: $MAX_KAKS (filters positive selection outliers)
  - Minimum alignment length: $MIN_ALIGNMENT_LENGTH codons

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i output/glycine_max/ks_results/ \\
     -o output/glycine_max/ks_summary.tsv

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
if [ -n "$INPUT_DIR" ] || [ -n "$OUTPUT_FILE" ]; then
    AUTO_DETECT=false
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    SPECIES_DIRS=(output/*/ks_results)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        INPUT_DIR="${SPECIES_DIRS[0]}"
        SPECIES_NAME=$(basename $(dirname "$INPUT_DIR"))
        OUTPUT_FILE="output/${SPECIES_NAME}/ks_summary.tsv"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Input:  $INPUT_DIR"
        echo "   Output: $OUTPUT_FILE"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species Ks result directories found"
        echo "       Specify paths explicitly with -i and -o"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No Ks results found in output/*/"
        echo "       Run pipeline_2/5_calculate_ks.sh first"
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
echo "   INPUT DIR   : $INPUT_DIR"
echo "   OUTPUT FILE : $OUTPUT_FILE"
echo "-- Quality filters:"
echo "   Max Ks      : $MAX_KS"
echo "   Max Ka      : $MAX_KA"
echo "   Max Ka/Ks   : $MAX_KAKS"
echo "   Min length  : $MIN_ALIGNMENT_LENGTH codons"

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    echo "       Run pipeline_2/5_calculate_ks.sh first"
    exit 1
fi

# Create output directory
mkdir -p "$(dirname "$OUTPUT_FILE")"

echo "-- Consolidating PAML yn00 results..."

# Count yn00 output files
TOTAL_FILES=$(find "$INPUT_DIR" -name "*_yn00.out" | wc -l)

echo "   Found $TOTAL_FILES yn00 output files"

if [ $TOTAL_FILES -eq 0 ]; then
    echo "ERROR: No yn00 output files (*_yn00.out) found in $INPUT_DIR"
    echo "       Run pipeline_2/5_calculate_ks.sh first"
    exit 1
fi

START_TIME=$(date +%s)

# Create header
echo -e "gene1\tgene2\tfamily\tks\tka\tka_ks\tlength\tstatus" > "$OUTPUT_FILE"

# Process all yn00 files and parse results
find "$INPUT_DIR" -name "*_yn00.out" | awk -v max_ks="$MAX_KS" -v max_ka="$MAX_KA" -v max_kaks="$MAX_KAKS" -v min_len="$MIN_ALIGNMENT_LENGTH" '
{
    file = $0
    
    # Extract gene names and family from path
    # Format: output/species/ks_results/family123/gene1_gene2_yn00.out
    split(file, path_parts, "/")
    family_dir = path_parts[length(path_parts)-1]  # family123
    filename = path_parts[length(path_parts)]      # gene1_gene2_yn00.out
    
    # Extract family number
    gsub(/family/, "", family_dir)
    family = family_dir
    
    # Extract gene names from filename
    gsub(/_yn00\.out$/, "", filename)
    split(filename, genes, "_")
    gene1 = genes[1]
    gene2 = genes[2]
    
    # Parse yn00 output file
    ks = "NA"
    ka = "NA"
    ka_ks = "NA"
    length_val = "NA"
    status = "FAILED"
    
    while ((getline line < file) > 0) {
        # Look for Yang & Nielsen 2000 method results
        if (line ~ /seq\. seq\.     S       N        t   kappa   omega/) {
            # This is the header line, skip empty line, get data line
            getline line < file  # skip empty line
            getline line < file  # get data line
            
            if (line != "" && line !~ /^$/) {
                # Format: seq1 seq2 S N t kappa omega dN +- SE dS +- SE
                split(line, fields)
                if (length(fields) >= 11) {
                    ka_ks = fields[7]   # omega (dN/dS)
                    ka = fields[8]      # dN value
                    ks = fields[11]     # dS value
                    status = "PARSED"
                }
            }
            break
        }
        
        # Get sequence length from ns/ls line
        if (line ~ /ns =.*ls =/) {
            # Extract ls value (sequence length in codons)
            if (match(line, /ls = *([0-9]+)/, arr)) {
                length_val = arr[1]
            }
        }
    }
    close(file)
    
    # Apply quality filters
    if (status == "PARSED") {
        # Check if values are numeric and within bounds
        if (ks != "NA" && ka != "NA" && ka_ks != "NA") {
            if (ks > max_ks) status = "KS_TOO_HIGH"
            else if (ka > max_ka) status = "KA_TOO_HIGH"  
            else if (ka_ks > max_kaks) status = "KAKS_TOO_HIGH"
            else if (length_val != "NA" && length_val < min_len) status = "TOO_SHORT"
            else status = "PASS"
        }
    }
    
    # Output result
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", gene1, gene2, family, ks, ka, ka_ks, length_val, status
    
    total++
    if (status == "PASS") passed++
    else if (status == "PARSED") filtered++
    else failed++
}

END {
    printf "-- Parsing summary:\n" > "/dev/stderr"
    printf "   Total files: %d\n", total > "/dev/stderr"
    printf "   Passed QC: %d\n", passed > "/dev/stderr"
    printf "   Filtered: %d\n", filtered > "/dev/stderr"  
    printf "   Failed: %d\n", failed > "/dev/stderr"
}
' >> "$OUTPUT_FILE"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

# Count results
TOTAL_PAIRS=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
PASSED_QC=$(tail -n +2 "$OUTPUT_FILE" | awk '$8=="PASS"' | wc -l)
FILTERED=$(tail -n +2 "$OUTPUT_FILE" | awk '$8!="PASS" && $8!="FAILED"' | wc -l)

echo ""
echo "-- Results consolidation completed"
echo "   Total pairs: $TOTAL_PAIRS"
echo "   Passed QC: $PASSED_QC"
echo "   Filtered: $FILTERED"
echo "   Processing time: ${ELAPSED}s"
echo "   Output file: $OUTPUT_FILE"

# Generate summary statistics
if [ $PASSED_QC -gt 0 ]; then
    echo ""
    echo "-- Summary statistics (QC-passed pairs):"
    
    tail -n +2 "$OUTPUT_FILE" | awk '$8=="PASS"' | awk -F'\t' '
    {
        if (NR==1) {
            min_ks = max_ks = $4
            min_ka = max_ka = $5
            min_kaks = max_kaks = $6
        }
        
        ks_sum += $4; ka_sum += $5; kaks_sum += $6
        
        if ($4 < min_ks) min_ks = $4
        if ($4 > max_ks) max_ks = $4
        if ($5 < min_ka) min_ka = $5
        if ($5 > max_ka) max_ka = $5
        if ($6 < min_kaks) min_kaks = $6
        if ($6 > max_kaks) max_kaks = $6
        
        count++
    }
    END {
        printf "   Ks: %.3f ± %.3f (range: %.3f - %.3f)\n", ks_sum/count, sqrt(ks_sum2/count - (ks_sum/count)^2), min_ks, max_ks
        printf "   Ka: %.3f ± %.3f (range: %.3f - %.3f)\n", ka_sum/count, sqrt(ka_sum2/count - (ka_sum/count)^2), min_ka, max_ka  
        printf "   Ka/Ks: %.3f ± %.3f (range: %.3f - %.3f)\n", kaks_sum/count, sqrt(kaks_sum2/count - (kaks_sum/count)^2), min_kaks, max_kaks
    }'
fi

echo ""
echo "-- Pipeline 2 completed successfully!"
echo "   Final results: $OUTPUT_FILE"

exit 0