#!/bin/bash

# -- keep only the longest isoform from the initial fasta file
set -euo pipefail

FASTA_FILE='data/peptides.fa'
FEATURE_FILE='data/protein_info.csv'

cat <<EOF
-- this script filters the initial fasta file to keep only the longest isoform per gene
   and creates:
     1) a filtered fasta file (peptides_longest.fa)
     2) a filtered protein info file (protein_info_longest.csv)
     
EOF

# -- get arguments if provided any
while getopts "f:i:h" flag; do
    case "${flag}" in
        f) FASTA_FILE="${OPTARG}" ;;
        i) FEATURE_FILE="${OPTARG}" ;;
        h)
            echo "Usage: $0 [-f FASTA_FILE] [-i FEATURE_FILE]"
            exit 0
            ;;
        *)
            echo "Invalid option. Use -h for help."
            exit 1
            ;;
    esac
done

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
echo "-- outputs:"
echo "   FASTA   : $FASTA_FILTERED"
echo "   FEATURES: $FEATURE_FILTERED"

# cleanup
rm -f "$TMP_METADATA" "$TMP_IDS"

echo "-- filtering completed successfully."
exit 0

