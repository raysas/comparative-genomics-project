#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    keeps only the longest isoform per gene from the initial peptide fasta file
#    and creates:
#      1) a filtered fasta file (e.g., peptides_longest.fa)
#      2) a filtered protein info file (e.g., protein_info_longest.csv)
# -- Usage:
#    bash ./pipeline_1/1_filter_isoforms.sh [-f FASTA_FILE] [-i FEATURE_FILE]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/1_filter_isoforms.sh -f "data/peptides.fa" -i "data/protein_info.csv"
# --------------------------------------------------------------------
###########################################################################

# ---------------------------------------------------------------------
# prepare variables, get arguments, set up logging
# ---------------------------------------------------------------------

# -- keep only the longest isoform from the initial fasta file
set -euo pipefail

# Default paths (will auto-detect species subdirectory)
FASTA_FILE=''
FEATURE_FILE=''
AUTO_DETECT=true

cat <<EOF
-- this script filters the initial fasta file to keep only the longest isoform per gene
   and creates:
     1) a filtered fasta file (peptides_longest.fa)
     2) a filtered protein info file (protein_info_longest.csv)
     
EOF

# -- get arguments if provided any
while getopts "f:i:h" flag; do
    case "${flag}" in
        f) FASTA_FILE="${OPTARG}"; AUTO_DETECT=false ;;
        i) FEATURE_FILE="${OPTARG}"; AUTO_DETECT=false ;;
        h)
            cat <<EOF
Usage: $0 [-f FASTA_FILE] [-i FEATURE_FILE]

OPTIONS:
  -f FASTA_FILE     Input peptide FASTA file
  -i FEATURE_FILE   Input protein info CSV file
  -h                Show this help message

AUTO-DETECTION:
  If no files are specified, the script will automatically detect
  species-specific directories under data/
  
  Expected structure:
    data/{species}/peptides.fa
    data/{species}/protein_info.csv

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -f data/glycine_max/peptides.fa -i data/glycine_max/protein_info.csv

EOF
            exit 0
            ;;
        *)
            echo "Invalid option. Use -h for help."
            exit 1
            ;;
    esac
done

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    # Look for data/species_name/ structure
    SPECIES_DIRS=(data/*/)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        SPECIES_DIR="${SPECIES_DIRS[0]%/}"
        SPECIES_NAME=$(basename "$SPECIES_DIR")
        
        FASTA_FILE="${SPECIES_DIR}/peptides.fa"
        FEATURE_FILE="${SPECIES_DIR}/protein_info.csv"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Using directory: $SPECIES_DIR"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species directories found in data/"
        echo "       Please specify files explicitly with -f and -i"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No species directories found in data/"
        echo "       Run 0_extract_data.sh first or specify files with -f and -i"
        exit 1
    fi
fi

# Validate files exist
if [ ! -f "$FASTA_FILE" ]; then
    echo "ERROR: FASTA file not found: $FASTA_FILE"
    exit 1
fi

if [ ! -f "$FEATURE_FILE" ]; then
    echo "ERROR: Feature file not found: $FEATURE_FILE"
    exit 1
fi

LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"


echo "-- filtering to keep only longest isoform per gene from:"
echo "   FASTA   : $FASTA_FILE"
echo "   FEATURES: $FEATURE_FILE"

# -- define output files
FASTA_DIR=$(dirname "$FASTA_FILE")
FEATURE_DIR=$(dirname "$FEATURE_FILE")

FASTA_FILTERED="${FASTA_DIR}/peptides_longest.fa"
FEATURE_FILTERED="${FEATURE_DIR}/protein_info_longest.csv"

TMP_METADATA="${FEATURE_DIR}/_tmp_longest_isoforms.csv"
TMP_IDS="${FEATURE_DIR}/_tmp_longest_ids.txt"

#########################################################################

# --------------------------------------------------------------------
# 1) select longest isoform per gene using FEATURE_FILE
# columns in protein_info.csv:
# peptide_id,gene_id,transcript_id,genome,chromosome,start_pos,end_pos,strand,description,length
# --------------------------------------------------------------------
echo "-- step 1: selecting longest isoform per gene based on 'length' column"

# header
head -n 1 "$FEATURE_FILE" > "$FEATURE_FILTERED"

# select best row (max length) per gene
tail -n +2 "$FEATURE_FILE" | awk -F',' '
{
    gene = $2
    len  = $10 + 0
    if (!(gene in best_len) || len > best_len[gene]) {
        best_len[gene] = len
        best_row[gene] = $0
    }
}
END {
    for (g in best_row) {
        print best_row[g]
    }
}' > "$TMP_METADATA"

cat "$TMP_METADATA" >> "$FEATURE_FILTERED"

# extract peptide IDs of kept isoforms
cut -d',' -f1 "$TMP_METADATA" > "$TMP_IDS"

N_GENES=$(wc -l < "$TMP_IDS" || echo 0)
echo "-- step 1 done: kept $N_GENES genes (one isoform each)"

# --------------------------------------------------------------------
# 2) filter FASTA to keep only selected peptide IDs
# --------------------------------------------------------------------
echo "-- step 2: filtering FASTA to keep only selected peptide IDs"
echo "   output FASTA: $FASTA_FILTERED"

python3 - "$FASTA_FILE" "$TMP_IDS" "$FASTA_FILTERED" << 'EOF'
import sys

fasta_path, ids_path, out_path = sys.argv[1:]

with open(ids_path) as f:
    keep_ids = {line.strip() for line in f if line.strip()}

with open(fasta_path) as fin, open(out_path, "w") as out:
    write = False
    for line in fin:
        if line.startswith(">"):
            pid = line[1:].split()[0]
            write = pid in keep_ids
        if write:
            out.write(line)
EOF

echo "-- step 2 done."
echo ""
echo "========================================="
echo "FILTERING COMPLETED"
echo "========================================="
echo "Filtered FASTA:    $FASTA_FILTERED"
echo "Filtered features: $FEATURE_FILTERED"
echo "Genes retained:    $N_GENES"
echo "========================================="

# cleanup
rm -f "$TMP_METADATA" "$TMP_IDS"

exit 0

