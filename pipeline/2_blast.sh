#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    1) creates a BLAST database from a given peptide fasta file
#    2) runs a BLASTP search of the peptides against themselves
# !! if the output directory already exists, it will skip running blast unless -f flag is provided
# -- Usage:
#    bash ./pipeline_1/2_blast.sh [-i INPUT_FILE] [-o OUTPUT_DIR] [-f] [-h]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/2_blast.sh -i "data/peptides_longest.fa" -o "output/blast_output"
# --------------------------------------------------------------------
###########################################################################

# ---------------------------------------------------------------------
# prepare variables, get arguments, set up logging
# ---------------------------------------------------------------------
# -- message on what this script does
cat <<EOF
-- this script:
  1) creates a BLAST database from a given peptide fasta file
  2) runs a BLASTP search of the peptides against themselves
  3) computes coverage and adds it to the last columns of the BLAST output
EOF

# -- default parameters
INPUT_FILE='data/peptides_longest.fa'
OUTPUT_DIR='output/blast_output'

DB_NAME='peptide_db'
OUTPUT_FILE='blast_results.tsv'
COV_OUTPUT_FILE="blast_results_with_coverage.tsv"


# -- make a flag to not keep blast_output directory across runs and make a new one each time
force_new=false

# -- get arguments
while getopts "i:o:f:h" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        f) force_new=true ;;
        h)
            echo "Usage: $0 [-i input_file] [-o output_directory] [-f] [-h]"
            echo "  -i    Input file"
            echo "  -o    Output directory"
            echo "  -f    Force creation of a new blast_output directory (do not reuse existing one - wont delete previous one will only rename it with prefix old_ and a number suffix)"
            echo "  -h    Show this help message"
            exit 0
            ;;
        *)
            echo "Invalid option: -${OPTARG}, follow the order of -i, -o, -f, -h just in case" >&2
            exit 1
            ;;
    esac
done


# -- name log file based on script name
LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo "-- We have the following parameters:"
echo "   INPUT FILE : $INPUT_FILE"
echo "   OUTPUT DIR : $OUTPUT_DIR"
echo "   FORCE NEW  : $force_new (create new output dir even if one exists)"


# -- check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist."
    exit 1
fi

# -------------------------------------------------------------------------
# -- main script logic
# -------------------------------------------------------------------------

# -- handle existing output directory
if [ "$force_new" = false ] && [ -d "$OUTPUT_DIR" ]; then
    echo "-- output directory $OUTPUT_DIR already exists and -f flag not set; using existing directory"
    echo "-- to force creation of a new output directory, use the -f flag"
# -- else output directory exists and -f flag is set: rename exsisting output dire to old_blast_output_1, old_blast_output_2, etc.
else
    echo "-- preparing to create new output directory: $OUTPUT_DIR"
    leaf_dir=$(basename "$OUTPUT_DIR")
    parent_dir=$(dirname "$OUTPUT_DIR")
    if [ -d "$OUTPUT_DIR" ]; then
        # -- move leaf to leaf_old_numbered
        n=1
        while [ -d "${parent_dir}/${leaf_dir}_old_${n}" ]; do
            ((n++))
        done
        mv "$OUTPUT_DIR" "${parent_dir}/${leaf_dir}_old_${n}"
        echo "-- ! moved existing directory $OUTPUT_DIR to ${parent_dir}/${leaf_dir}_old_${n}"
    fi

    # -- create output directory
    mkdir -p "$OUTPUT_DIR"
    echo "-- created output directory: $OUTPUT_DIR"

    # -- make blast db
    if ! command -v makeblastdb &> /dev/null
    then
        echo "-- error: makeblastdb could not be found: install BLAST+ tools"
        exit 2
    fi

    makeblastdb -in "$INPUT_FILE" -dbtype prot -out "${OUTPUT_DIR}/${DB_NAME}"
    echo "-- created BLAST database from $INPUT_FILE"

    # -- run blastp
    if ! command -v blastp &> /dev/null
    then
        echo "-- error: blastp could not be found: install BLAST+ tools"
        exit 3
    fi

    blastp -query "$INPUT_FILE" -db "${OUTPUT_DIR}/${DB_NAME}" -out "${OUTPUT_DIR}/${OUTPUT_FILE}" -outfmt 6 -max_target_seqs 1
    echo "-- BLASTP search completed, results saved to ${OUTPUT_DIR}/${OUTPUT_FILE}"
    current_timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "-- completed at: $current_timestamp"
fi


exit 0
