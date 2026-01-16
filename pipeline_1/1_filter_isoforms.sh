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

<<<<<<< HEAD
# Default paths (will auto-detect species subdirectory)
FASTA_FILE=''
FEATURE_FILE=''
AUTO_DETECT=true
=======
# Resolve paths to run from anywhere and improve UX
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
>>>>>>> master

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
NC='\033[0m'

# Default paths (can auto-detect species subdirectory)
FASTA_FILE=''
FEATURE_FILE=''
AUTO_DETECT=false

# cat <<EOF
# -- this script filters the initial fasta file to keep only the longest isoform per gene
#    and creates:
#      1) a filtered fasta file (peptides_longest.fa)
#      2) a filtered protein info file (protein_info_longest.csv)
     
# EOF

# -- get arguments if provided any
while getopts "f:i:s:h" flag; do
    case "${flag}" in
<<<<<<< HEAD
        f) FASTA_FILE="${OPTARG}"; AUTO_DETECT=false ;;
        i) FEATURE_FILE="${OPTARG}"; AUTO_DETECT=false ;;
        h)
            cat <<EOF
Usage: $0 [-f FASTA_FILE] [-i FEATURE_FILE]
=======
        f) FASTA_FILE="${OPTARG}";;
        i) FEATURE_FILE="${OPTARG}";;
        s) SPECIES_NAME="${OPTARG}"; AUTO_DETECT=true ;;
        h)
            cat <<EOF
Usage: $0 [-f FASTA_FILE] [-i FEATURE_FILE] [-s SPECIES_NAME] [-h]
>>>>>>> master

OPTIONS:
  -f FASTA_FILE     Input peptide FASTA file
  -i FEATURE_FILE   Input protein info CSV file
<<<<<<< HEAD
=======
  -s SPECIES_NAME   Specify species name (for auto-detection)
>>>>>>> master
  -h                Show this help message

AUTO-DETECTION:
  If no files are specified, the script will automatically detect
  species-specific directories under data/
  
<<<<<<< HEAD
  Expected structure:
    data/{species}/peptides.fa
    data/{species}/protein_info.csv

EXAMPLES:
  # Auto-detect (recommended)
  $0
=======
    Expected structure (relative to repo root):
        data/{species}/peptides.fa
        data/{species}/protein_info.csv

EXAMPLES:
  # Auto-detect (recommended)
  $0 -s glycine_max
>>>>>>> master

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
<<<<<<< HEAD
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
=======
    echo -e "${YELLOW}-- Auto-detecting species directory...${NC}"
    if [ -n "${SPECIES_NAME:-}" ]; then
        SPECIES_DIR="$REPO_ROOT/data/${SPECIES_NAME}"
        if [ -d "$SPECIES_DIR" ]; then
            FASTA_FILE="${SPECIES_DIR}/peptides.fa"
            FEATURE_FILE="${SPECIES_DIR}/processed/protein_info.csv"
            echo -e "   Detected species: ${GREEN}$SPECIES_NAME${NC}"
            echo -e "   Using directory: ${BLUE}$SPECIES_DIR${NC}"
        else
            echo -e "${RED}ERROR: Specified species directory not found: $SPECIES_DIR${NC}"
            exit 1
        fi
    fi
else
    echo -e "${RED}ERROR: Species name not specified. Use -s <species> or provide -f and -i.${NC}"
    exit 1
>>>>>>> master
fi

# Fallback defaults if not provided (run-from-anywhere)
if [ -z "$FASTA_FILE" ] && [ -z "$FEATURE_FILE" ]; then
    if [ -f "$REPO_ROOT/data/peptides.fa" ] && [ -f "$REPO_ROOT/data/protein_info.csv" ]; then
        FASTA_FILE="$REPO_ROOT/data/peptides.fa"
        FEATURE_FILE="$REPO_ROOT/data/protein_info.csv"
        echo -e "${YELLOW}-- Using defaults from repo root: data/peptides.fa, data/protein_info.csv${NC}"
    fi
fi

# Normalize to absolute paths if possible
if command -v realpath >/dev/null 2>&1; then
    [ -n "$FASTA_FILE" ] && FASTA_FILE="$(realpath "$FASTA_FILE")"
    [ -n "$FEATURE_FILE" ] && FEATURE_FILE="$(realpath "$FEATURE_FILE")"
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

# Setup logging (always inside the pipeline directory)
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh)_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN}  Step 1: LONGEST ISOFORM FILTER ${NC}"
echo -e "${GREEN}=========================================${NC}"
echo -e "Input Files:"
echo -e "  FASTA    : ${BLUE}$FASTA_FILE${NC}"
echo -e "  Features : ${BLUE}$FEATURE_FILE${NC}"

# -- define output files
FASTA_DIR=$(dirname "$FASTA_FILE")
FEATURE_DIR=$(dirname "$FEATURE_FILE")

FASTA_FILTERED="${FEATURE_DIR}/peptides_longest.fa"
FEATURE_FILTERED="${FEATURE_DIR}/protein_info_longest.csv"

TMP_METADATA="${FEATURE_DIR}/_tmp_longest_isoforms.csv"
TMP_IDS="${FEATURE_DIR}/_tmp_longest_ids.txt"

#########################################################################

# --------------------------------------------------------------------
# 1) select longest isoform per gene using FEATURE_FILE
# columns in protein_info.csv:
# peptide_id,gene_id,transcript_id,genome,chromosome,start_pos,end_pos,strand,description,length
# --------------------------------------------------------------------
# echo -e "${YELLOW}Step 1:${NC} Selecting longest isoform per gene (by length)"

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
# echo -e "${GREEN}✓ Step 1 complete:${NC} kept $N_GENES genes (one isoform each)"

# --------------------------------------------------------------------
# 2) filter FASTA to keep only selected peptide IDs
# --------------------------------------------------------------------
# echo -e "${YELLOW}Step 2:${NC} Filtering FASTA to keep selected peptide IDs"
echo -e "  Output FASTA: ${BLUE}$FASTA_FILTERED${NC}"

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

<<<<<<< HEAD
echo "-- step 2 done."
echo ""
echo "========================================="
echo "FILTERING COMPLETED"
echo "========================================="
echo "Filtered FASTA:    $FASTA_FILTERED"
echo "Filtered features: $FEATURE_FILTERED"
echo "Genes retained:    $N_GENES"
echo "========================================="
=======
# echo -e "${GREEN}✓ Step 2 complete${NC}"
# echo ""
echo -e "${GREEN}=========================================${NC}"
echo -e "${GREEN} FILTERING COMPLETED ${NC}"
echo -e "${GREEN}=========================================${NC}"
echo -e "Filtered FASTA   : ${BLUE}$FASTA_FILTERED${NC}"
echo -e "Filtered features: ${BLUE}$FEATURE_FILTERED${NC}"
echo -e "Genes retained   : ${YELLOW}$N_GENES${NC}"
echo -e "Log              : ${BLUE}$LOG_FILE${NC}"
echo -e "${GREEN}=========================================${NC}"
>>>>>>> master

# cleanup
rm -f "$TMP_METADATA" "$TMP_IDS"

exit 0

