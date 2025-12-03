#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    Takes gene family clusters from pipeline_1 and creates pairwise 
#    protein FASTA files for each family (needed for Ks calculation)
# -- Usage:
#    bash ./pipeline_2/1_prepare_pairs.sh [-i FAMILY_FILE] [-p PROTEIN_FASTA] [-o OUTPUT_DIR] [-h]
# -- default (without params) equivalent to:
#    Auto-detects from species directories
# --------------------------------------------------------------------

# -- default parameters
FAMILY_FILE=''
PROTEIN_FASTA=''
OUTPUT_DIR=''
AUTO_DETECT=true
NUM_JOBS=$(nproc)

# -- arguments
while getopts "i:p:o:h" flag; do
    case "${flag}" in
        i) FAMILY_FILE="${OPTARG}" ;;
        p) PROTEIN_FASTA="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE       Input family clustering file (.tsv)
  -p FILE       Protein FASTA file (peptides_longest.fa)
  -o DIR        Output directory for family pairs
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
    Input families: output/{species}/clusters/protein_families_*.tsv
    Protein FASTA:  data/{species}/peptides_longest.fa
    Output:         output/{species}/families/

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i output/glycine_max/clusters/protein_families_network_50families_max5.tsv \\
     -p data/glycine_max/peptides_longest.fa \\
     -o output/glycine_max/families/

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
if [ -n "$FAMILY_FILE" ] || [ -n "$PROTEIN_FASTA" ] || [ -n "$OUTPUT_DIR" ]; then
    AUTO_DETECT=false
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    SPECIES_DIRS=(output/*/clusters)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        CLUSTER_DIR="${SPECIES_DIRS[0]}"
        SPECIES_NAME=$(basename $(dirname "$CLUSTER_DIR"))
        
        # Find family file (prefer test dataset if available)
        if [ -f "${CLUSTER_DIR}/protein_families_test_dataset.tsv" ]; then
            FAMILY_FILE="${CLUSTER_DIR}/protein_families_test_dataset.tsv"
            echo "   Using balanced test dataset"
        elif [ -f "${CLUSTER_DIR}/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv" ]; then
            FAMILY_FILE="${CLUSTER_DIR}/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv"
            echo "   Using full dataset"
        else
            echo "ERROR: No family clustering file found in $CLUSTER_DIR"
            echo "       Expected files:"
            echo "         - ${CLUSTER_DIR}/protein_families_test_dataset.tsv (for testing)"
            echo "         - ${CLUSTER_DIR}/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv (full dataset)"
            echo "       Run pipeline_1 clustering or create test dataset with:"
            echo "         ./create_test_dataset.sh"
            exit 1
        fi
        
        PROTEIN_FASTA="data/${SPECIES_NAME}/peptides_longest.fa"
        OUTPUT_DIR="output/${SPECIES_NAME}/families"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Family file: $FAMILY_FILE"
        echo "   Protein FASTA: $PROTEIN_FASTA"
        echo "   Output: $OUTPUT_DIR"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species cluster directories found"
        echo "       Specify files explicitly with -i, -p, -o"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No cluster results found in output/*/"
        echo "       Run pipeline_1 (steps 1-6) first"
        exit 1
    fi
fi

# Setup logging
LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then

    PAIR_LIST="/tmp/pair_list.tsv"
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo "-- Parameters:"
echo "   FAMILY FILE    : $FAMILY_FILE"
echo "   PROTEIN FASTA  : $PROTEIN_FASTA"
echo "   OUTPUT DIR     : $OUTPUT_DIR"

# Check input files exist
if [ ! -f "$FAMILY_FILE" ]; then
    echo "ERROR: Family file not found: $FAMILY_FILE"
    echo "       Run pipeline_1 clustering step first"
    exit 1
fi

if [ ! -f "$PROTEIN_FASTA" ]; then
    echo "ERROR: Protein FASTA file not found: $PROTEIN_FASTA"
    echo "       Run pipeline_1 steps 0-1 first"
    exit 1
fi

echo "-- Preparing pairwise protein files from gene families..."

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count families for progress
TOTAL_FAMILIES=$(tail -n +2 "$FAMILY_FILE" | cut -f2 | sort -u | wc -l)
echo "   Processing $TOTAL_FAMILIES families..."

START_TIME=$(date +%s)

echo "-- Pre-creating family directories..."
# tail -n +2 "$FAMILY_FILE" | cut -f2 | sort -u | while read fam; do
#     clean_fam=$(echo "$fam" | tr -d '\r' | tr -d '[:space:]')
#     mkdir -p "$OUTPUT_DIR/family${clean_fam}"
# done

# Pre-create all family directories
# find "/tmp/pipeline2_full/alignments/" -maxdepth 1 -type d -name "family*" -printf "%f\n" | \
    # parallel -j "$NUM_JOBS" "mkdir -p '$OUTPUT_DIR/{}'"


echo "-- Generating the list of pairwise FASTA files..."
PAIR_LIST="/tmp/pair_list.tsv"
# awk -F'\t' -v protein_file="$PROTEIN_FASTA" -v output_dir="$OUTPUT_DIR" '
# BEGIN {
#     while ((getline line < protein_file) > 0) {
#         if (line ~ /^>/) {
#             split(substr(line, 2), parts, " ")
#             seq_id = parts[1]
#         } else {
#             sequences[seq_id] = sequences[seq_id] line
#         }
#     }
#     close(protein_file)
# }
# NR == 1 { next }
# {
#     gene_id = $1
#     family_id = $2
#     gsub(/\r/, "", gene_id)
#     gsub(/\r/, "", family_id)
#     gsub(/ /, "", gene_id)
#     gsub(/ /, "", family_id)
#     family_genes[family_id][++family_count[family_id]] = gene_id
# }
# END {
#     for (fam in family_genes) {
#         clean_fam = fam
#         gsub(/\r/, "", clean_fam)
#         gsub(/ /, "", clean_fam)
#         fam_dir = output_dir "family" clean_fam
#         for (i = 1; i <= family_count[fam]; i++) {
#             for (j = i + 1; j <= family_count[fam]; j++) {
#                 gene1 = family_genes[fam][i]
#                 gene2 = family_genes[fam][j]
#                 if (gene1 in sequences && gene2 in sequences) {
#                     print fam_dir "\t" gene1 "\t" gene2 "\t" sequences[gene1] "\t" sequences[gene2]
#                 }
#             }
#         }
#     }
# }
# ' "$FAMILY_FILE" > "$PAIR_LIST"

# Step 2: GNU Parallel to create FASTA files with progress bar
# TOTAL_PAIRS=$(wc -l < "$PAIR_LIST")

# echo "Creating $TOTAL_PAIRS pairwise FASTA files in parallel..."

function make_fasta {
    fam_dir="$1"; gene1="$2"; gene2="$3"; seq1="$4"; seq2="$5"
    out_file="$fam_dir/${gene1}_${gene2}.fa"
    if [ -f "$out_file" ] && [ -s "$out_file" ]; then
        return 0
    fi
    echo ">$gene1" > "$out_file"
    echo "$seq1" >> "$out_file"
    echo ">$gene2" >> "$out_file"
    echo "$seq2" >> "$out_file"
}
export -f make_fasta


# Filter out lines that would exceed ARG_MAX for command line length
ARG_MAX=25000
tail -n +430000 "$PAIR_LIST" | awk -F'\t' -v max="$ARG_MAX" '{sum=0; for(i=1;i<=NF;i++) sum+=length($i); if(sum+NF-1 < max) print $0}' | \
    parallel -j "$NUM_JOBS" --colsep '\t' --progress make_fasta {1} {2} {3} {4} {5}

# rm -f "$PAIR_LIST"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "-- Pair preparation completed successfully"
echo "   Processing time: ${ELAPSED}s"
echo "   Output directory: $OUTPUT_DIR"
echo "   Family subdirectories: output/{species}/pairs/family*/"

echo "-- Next step: Run protein alignment on families"
echo "   bash pipeline_2/3_align_proteins.sh"

exit 0