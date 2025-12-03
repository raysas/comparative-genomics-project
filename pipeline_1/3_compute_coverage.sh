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
cat <<EOF
-- this script computes query and subject coverage from BLASTP output
   and appends the coverage values as new columns to the BLAST output file
EOF

# -- default parameters
INPUT_FILE=''
OUTPUT_FILE=''
PROTEIN_INFO_FILE=''
AUTO_DETECT=true

# -- arguments
while getopts "i:o:p:fh" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_FILE="${OPTARG}" ;;
        p) PROTEIN_INFO_FILE="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE       Input BLAST results file (blast_results.tsv)
  -o FILE       Output file with coverage (blast_results_with_coverage.tsv)
  -p FILE       Protein info file (protein_info_longest.csv)
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
    Input:  output/{species}/blast_output/blast_results.tsv
    Output: output/{species}/blast_output/blast_results_with_coverage.tsv
    Protein: data/{species}/protein_info_longest.csv

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
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo "-- Parameters:"
echo "   INPUT FILE    : $INPUT_FILE"
echo "   PROTEIN INFO  : $PROTEIN_INFO_FILE"
echo "   OUTPUT FILE   : $OUTPUT_FILE"

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

echo "-- Computing coverage for BLASTP results..."

################################################################################
# Optimized coverage computation:
# Single-pass AWK script that loads protein lengths into memory hash table
# and processes BLAST results in one streaming operation (no temp files)
################################################################################

# Count input rows for progress
TOTAL_BLAST_HITS=$(wc -l < "$INPUT_FILE")
echo "   Processing $TOTAL_BLAST_HITS BLAST hits..."

# Start timing
START_TIME=$(date +%s)

# Single optimized AWK script - no temp files, no joins, no sorts
echo "   Computing coverage with optimized single-pass algorithm..."

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
    printf "   Loaded %d protein lengths into memory\n", length(protein_lengths) > "/dev/stderr"
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
    
    # Skip if either protein not found
    if (qlength == "" || slength == "") {
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

echo ""
echo "-- Coverage computation completed successfully"
echo "   Input hits:    $TOTAL_BLAST_HITS"
echo "   Output hits:   $FINAL_HITS" 
echo "   Processing time: ${ELAPSED}s"
echo "   Output file:   $OUTPUT_FILE"

# Validation
if [ $FINAL_HITS -ne $TOTAL_BLAST_HITS ]; then
    MISSING=$((TOTAL_BLAST_HITS - FINAL_HITS))
    echo "   WARNING: $MISSING hits missing (proteins not in protein_info file)"
fi

exit 0