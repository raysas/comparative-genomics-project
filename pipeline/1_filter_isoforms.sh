#!/bin/bash

# -- keep only the longest isoform from the initial fasta file

FASTA_FILE='data/peptides.fa'
FEATURE_FILE='data/protein_info.csv'

LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1

cat <<EOF
-- this script filters the initial fasta file to keep only the longest isoform per gene
   and updates the protein info file accordingly
EOF

# -- get arguments if provided any
while getopts f:i:h flag
do
    case "${flag}" in
        f) FASTA_FILE=${OPTARG};;
        i) FEATURE_FILE=${OPTARG};;
        h) echo "Usage: $0 [-f FASTA_FILE] [-i FEATURE_FILE]"
           exit 0;;
    esac
done


echo "-- filtering to keep only longest isoform per gene from $FASTA_FILE"

