#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    computes query and subject coverage from BLASTP output
#    and appends the coverage values as new columns to the BLAST output file
#
# -- splitted into 2 steps:
# I) preparing BLAST output with query and subject lengths
# II) computing coverage and appending to BLAST output
# -- Usage:
#    bash ./pipeline_1/3_compute_coverage.sh [-i INPUT_FILE] [-o COVERAGE_OUTPUT_FILE] [-p PROTEIN_INFO_FILE] [-h]
# -- default (without params) equivalent to:
#    Auto-detects from species directories
# --------------------------------------------------------------------
###########################################################################

# ---------------------------------------------------------------------
# prepare variables, get arguments, set up logging
# ---------------------------------------------------------------------
# -- message on what this script does
# Resolve paths to run from anywhere
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
NC='\033[0m'

echo -e "${GREEN}====================================${NC}"
echo -e "${GREEN} Step 3: Coverage ${NC}"
echo -e "${GREEN}====================================${NC}"
# echo -e " 1) Prepare BLAST output with lengths"
# echo -e " 2) Compute query/subject coverage and append"

# -- default parameters
INPUT_FILE=''
OUTPUT_FILE=''
PROTEIN_INFO_FILE=''
AUTO_DETECT=true
<<<<<<< HEAD

# -- arguments
while getopts "i:o:p:fh" flag; do
=======
SPECIES_NAME=""

# -- arguments
while getopts "i:o:p:s:h" flag; do
>>>>>>> master
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_FILE="${OPTARG}" ;;
        p) PROTEIN_INFO_FILE="${OPTARG}" ;;
        s) SPECIES_NAME="${OPTARG}"; AUTO_DETECT=true ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE       Input BLAST results file (blast_results.tsv)
  -o FILE       Output file with coverage (blast_results_with_coverage.tsv)
  -p FILE       Protein info file (protein_info_longest.csv)
<<<<<<< HEAD
=======
    -s NAME       Species name (uses data/NAME and output/NAME paths)
>>>>>>> master
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
<<<<<<< HEAD
    Input:  output/{species}/blast_output/blast_results.tsv
    Output: output/{species}/blast_output/blast_results_with_coverage.tsv
    Protein: data/{species}/protein_info_longest.csv
=======
        Input:  output/pipeline1/{species}/blast_results/blast_results.tsv
        Output: output/pipeline1/{species}/blast_results/blast_results_with_coverage.tsv
        Protein: data/{species}/protein_info_longest.csv
>>>>>>> master

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i output/glycine_max/blast_output/blast_results.tsv \\
     -o output/glycine_max/blast_output/blast_results_with_coverage.tsv \\
     -p data/glycine_max/protein_info_longest.csv

EOF
            exit 0
            ;;
        *)
            echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
    esac
done

# Check if user provided any input - disable auto-detect only if files are specified
if [ -n "$INPUT_FILE" ] || [ -n "$OUTPUT_FILE" ] || [ -n "$PROTEIN_INFO_FILE" ]; then
    AUTO_DETECT=false
<<<<<<< HEAD
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    SPECIES_DIRS=(output/*/blast_output)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        BLAST_DIR="${SPECIES_DIRS[0]}"
        SPECIES_NAME=$(basename $(dirname "$BLAST_DIR"))
        
        INPUT_FILE="${BLAST_DIR}/blast_results.tsv"
        OUTPUT_FILE="${BLAST_DIR}/blast_results_with_coverage.tsv"
        PROTEIN_INFO_FILE="data/${SPECIES_NAME}/protein_info_longest.csv"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Input:  $INPUT_FILE"
        echo "   Output: $OUTPUT_FILE"
        echo "   Protein info: $PROTEIN_INFO_FILE"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species BLAST results found"
        echo "       Specify files explicitly with -i, -o, -p"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No BLAST results found in output/*/"
        echo "       Run 2_blast.sh first"
        exit 1
    fi
fi

# Set default output file name if not provided or if OUTPUT_FILE is a directory
if [ -z "$OUTPUT_FILE" ] && [ -n "$INPUT_FILE" ]; then
    OUTPUT_DIR=$(dirname "$INPUT_FILE")
    OUTPUT_FILE="${OUTPUT_DIR}/$(basename "$INPUT_FILE" .tsv)_with_coverage.tsv"
elif [ -d "$OUTPUT_FILE" ]; then
    # If OUTPUT_FILE is a directory, append default filename
    OUTPUT_FILE="${OUTPUT_FILE%/}/$(basename "$INPUT_FILE" .tsv)_with_coverage.tsv"
fi

LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
=======
>>>>>>> master
fi


# Auto-detect species directory (aligned with Step 2 layout)
if [ "$AUTO_DETECT" = true ]; then
    echo -e "${YELLOW}-- Auto-detecting species directory...${NC}"
    if [ -n "$SPECIES_NAME" ]; then
        SPECIES_DIR="$REPO_ROOT/data/$SPECIES_NAME"
        BLAST_DIR="$REPO_ROOT/output/pipeline1/${SPECIES_NAME}/blast_results"
        if [ -d "$SPECIES_DIR" ] && [ -d "$BLAST_DIR" ]; then
            INPUT_FILE="${BLAST_DIR}/blast_results.tsv"
            OUTPUT_FILE="${BLAST_DIR}/blast_results_with_coverage.tsv"
            PROTEIN_INFO_FILE="${SPECIES_DIR}/processed/protein_info_longest.csv"
            [ ! -f "$PROTEIN_INFO_FILE" ] && PROTEIN_INFO_FILE="${SPECIES_DIR}/protein_info_longest.csv"
            echo -e "   Detected species: ${GREEN}$SPECIES_NAME${NC}"
            # echo -e "   Input : ${BLUE}$INPUT_FILE${NC}"
            # echo -e "   Output: ${BLUE}$OUTPUT_FILE${NC}"
            # echo -e "   Protein info: ${BLUE}$PROTEIN_INFO_FILE${NC}"
        else
            echo -e "${RED}ERROR:${NC} Species directories not found for '$SPECIES_NAME'"
            exit 1
        fi
    else
        SPECIES_BLAST_DIRS=("$REPO_ROOT"/output/pipeline1/*/blast_results)
        if [ ${#SPECIES_BLAST_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_BLAST_DIRS[0]}" ]; then
            BLAST_DIR="${SPECIES_BLAST_DIRS[0]}"
            SPECIES_NAME=$(basename $(dirname "$BLAST_DIR"))
            SPECIES_DIR="$REPO_ROOT/data/$SPECIES_NAME"
            INPUT_FILE="${BLAST_DIR}/blast_results.tsv"
            OUTPUT_FILE="${BLAST_DIR}/blast_results_with_coverage.tsv"
            PROTEIN_INFO_FILE="${SPECIES_DIR}/processed/protein_info_longest.csv"
            [ ! -f "$PROTEIN_INFO_FILE" ] && PROTEIN_INFO_FILE="${SPECIES_DIR}/protein_info_longest.csv"
            echo -e "   Detected species: ${GREEN}$SPECIES_NAME${NC}"
            echo -e "   Input : ${BLUE}$INPUT_FILE${NC}"
            echo -e "   Output: ${BLUE}$OUTPUT_FILE${NC}"
            echo -e "   Protein info: ${BLUE}$PROTEIN_INFO_FILE${NC}"
        elif [ ${#SPECIES_BLAST_DIRS[@]} -gt 1 ]; then
            echo -e "${RED}ERROR:${NC} Multiple species BLAST results found. Use -s or specify -i/-o/-p."
            exit 1
        else
            echo -e "${RED}ERROR:${NC} No BLAST results found. Run 2_blast.sh first."
            exit 1
        fi
    fi
fi

# Set default output file name if not provided or if OUTPUT_FILE is a directory
if [ -z "$OUTPUT_FILE" ] && [ -n "$INPUT_FILE" ]; then
    OUTPUT_DIR=$(dirname "$INPUT_FILE")
    OUTPUT_FILE="${OUTPUT_DIR}/$(basename "$INPUT_FILE" .tsv)_with_coverage.tsv"
elif [ -d "$OUTPUT_FILE" ]; then
    # If OUTPUT_FILE is a directory, append default filename
    OUTPUT_FILE="${OUTPUT_FILE%/}/$(basename "$INPUT_FILE" .tsv)_with_coverage.tsv"
fi

LOG_DIR="${SCRIPT_DIR}/logs/pipeline"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh)_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"
echo -e "${GREEN}Input:   ${BLUE}$INPUT_FILE${NC}"
echo -e "${GREEN}Protein: ${BLUE}$PROTEIN_INFO_FILE${NC}"
echo -e "${GREEN}Output:  ${BLUE}$OUTPUT_FILE${NC}"

<<<<<<< HEAD
echo "-- Parameters:"
echo "   INPUT FILE    : $INPUT_FILE"
echo "   PROTEIN INFO  : $PROTEIN_INFO_FILE"
echo "   OUTPUT FILE   : $OUTPUT_FILE"

=======
>>>>>>> master
# -- check input files exist
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: BLAST results file not found: $INPUT_FILE"
    echo "       Run 2_blast.sh first"
    exit 1
fi

if [ ! -f "$PROTEIN_INFO_FILE" ]; then
    echo "ERROR: Protein info file not found: $PROTEIN_INFO_FILE"
    echo "       Run 0_extract_data.sh and 1_filter_isoforms.sh first"
    exit 1
fi

# Check if output already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "-- Output file already exists: $OUTPUT_FILE"
    echo "   Overwriting existing file"
fi

<<<<<<< HEAD
echo "-- Computing coverage for BLASTP results..."
=======
echo -e "${YELLOW}-- Computing coverage for BLASTP results...${NC}"
>>>>>>> master

################################################################################
# Optimized coverage computation:
# Single-pass AWK script that loads protein lengths into memory hash table
# and processes BLAST results in one streaming operation (no temp files)
################################################################################

# Count input rows for progress
TOTAL_BLAST_HITS=$(wc -l < "$INPUT_FILE")
<<<<<<< HEAD
echo "   Processing $TOTAL_BLAST_HITS BLAST hits..."
=======
echo -e "   Processing: ${YELLOW}$TOTAL_BLAST_HITS${NC} BLAST hits"
>>>>>>> master

# Start timing
START_TIME=$(date +%s)

# Single optimized AWK script - no temp files, no joins, no sorts
<<<<<<< HEAD
echo "   Computing coverage with optimized single-pass algorithm..."
=======
echo -e "   Using optimized single-pass AWK"
>>>>>>> master

# Write header first
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlength\tslength\tqcov\tscov" > "$OUTPUT_FILE"

# Process everything in one AWK command
awk -F'\t' -v protein_file="$PROTEIN_INFO_FILE" 'BEGIN {
    OFS = "\t"
    line_count = 0
    processed = 0
    missing = 0
    
    # Load all protein lengths into memory hash table
    while ((getline line < protein_file) > 0) {
        line_count++
        if (line_count == 1) continue  # skip header
        split(line, fields, ",")
        if (length(fields) >= 10) {
            protein_id = fields[1]
            length_val = fields[10]
            protein_lengths[protein_id] = length_val
        }
    }
    close(protein_file)
<<<<<<< HEAD
    printf "   Loaded %d protein lengths into memory\n", length(protein_lengths) > "/dev/stderr"
=======
    printf "   Loaded %d protein lengths\n", length(protein_lengths) > "/dev/stderr"
>>>>>>> master
}
{
    # Process BLAST results
    qseqid = $1; sseqid = $2
    pident = $3; aln_length = $4; mismatch = $5; gapopen = $6
    qstart = $7; qend = $8; sstart = $9; send = $10
    evalue = $11; bitscore = $12
    
    # Look up lengths from hash table (O(1) operation)
    qlength = protein_lengths[qseqid]
    slength = protein_lengths[sseqid]
    
<<<<<<< HEAD
    # Skip if either protein not found
    if (qlength == "" || slength == "") {
=======
    # Skip if either protein not found or length is not a positive integer
    if (qlength == "" || slength == "" || qlength <= 0 || slength <= 0) {
>>>>>>> master
        missing++
        next
    }
    
    # Calculate coverage percentages
    qcov = (qend - qstart + 1) / qlength * 100
    scov = (send - sstart + 1) / slength * 100
    
    # Output result
    print qseqid, sseqid, pident, aln_length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlength, slength, qcov, scov
    processed++
}
END {
    printf "   Processed: %d hits\n", processed > "/dev/stderr"
    if (missing > 0) {
        printf "   WARNING: Skipped %d hits due to missing protein lengths\n", missing > "/dev/stderr"
    }
}' "$INPUT_FILE" >> "$OUTPUT_FILE"

# Calculate timing and final stats
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
FINAL_HITS=$(tail -n +2 "$OUTPUT_FILE" | wc -l)

<<<<<<< HEAD
echo ""
echo "-- Coverage computation completed successfully"
echo "   Input hits:    $TOTAL_BLAST_HITS"
echo "   Output hits:   $FINAL_HITS" 
echo "   Processing time: ${ELAPSED}s"
echo "   Output file:   $OUTPUT_FILE"

# Validation
if [ $FINAL_HITS -ne $TOTAL_BLAST_HITS ]; then
    MISSING=$((TOTAL_BLAST_HITS - FINAL_HITS))
=======
echo -e "\n${GREEN}✓ Coverage computation completed${NC}"
echo -e "   Input hits:    ${BLUE}$TOTAL_BLAST_HITS${NC}"
echo -e "   Output hits:   ${BLUE}$FINAL_HITS${NC}"
echo -e "   Time:          ${YELLOW}${ELAPSED}s${NC}"
echo -e "   Output file:   ${BLUE}$OUTPUT_FILE${NC}"


# Validation (robust arithmetic and quoting)
if [ -n "$FINAL_HITS" ] && [ -n "$TOTAL_BLAST_HITS" ] && [ "$FINAL_HITS" -ne "$TOTAL_BLAST_HITS" ]; then
    if [ "$FINAL_HITS" -gt "$TOTAL_BLAST_HITS" ]; then
        MISSING=0
    else
        MISSING=$((TOTAL_BLAST_HITS - FINAL_HITS))
    fi
>>>>>>> master
    echo "   WARNING: $MISSING hits missing (proteins not in protein_info file)"
fi

exit 0