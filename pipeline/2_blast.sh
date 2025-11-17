#!/bin/bash

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


# -- name log file based on script name
LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1


# -- get input arg from command line
while getopts i:o:h flag
do
    case "${flag}" in
        i) INPUT_FILE=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        h) echo "Usage: $0 [-i input_file] [-o output_directory]"
           exit 0;;
    esac
done

# -- check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file $INPUT_FILE does not exist."
    exit 1
fi

# -- handle existing output directory
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

blastp -query "$INPUT_FILE" -db "${OUTPUT_DIR}/${DB_NAME}" -out "${OUTPUT_DIR}/${OUTPUT_FILE}" -outfmt 6 
echo "-- BLASTP search completed, results saved to ${OUTPUT_DIR}/${OUTPUT_FILE}"

compute_cov() {
    awk '{
        q_start = $7
        q_end   = $8
        s_start = $9
        s_end   = $10
        q_length = $13
        s_length = $15
        id = $3
        bit = $12

        q_cov = (q_end - q_start + 1) / q_length
        s_cov = (s_end - s_start + 1) / s_length

        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, q_cov, s_cov
    }' "$1"
}
# -- compute coverage and add to BLAST output
compute_cov "${OUTPUT_DIR}/${OUTPUT_FILE}" > "${OUTPUT_DIR}/${COV_OUTPUT_FILE}"
echo "-- computed coverage and saved to ${OUTPUT_DIR}/${COV_OUTPUT_FILE}"
echo "   last 2 columns correspond to qcoverage and scoverage"

# -- fix header based on qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcov scov
echo "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcov	scov" > "${OUTPUT_DIR}/header.txt"
cat "${OUTPUT_DIR}/header.txt" "${OUTPUT_DIR}/${COV_OUTPUT_FILE}" > "${OUTPUT_DIR}/blast_results_final.tsv"
mv "${OUTPUT_DIR}/blast_results_final.tsv" "${OUTPUT_DIR}/${COV_OUTPUT_FILE}"
rm "${OUTPUT_DIR}/header.txt"


echo "-- BLAST pipeline finished successfully"

exit 0
