#!/bin/bash

# -- splitted into 2 steps:
# I) preparing BLAST output with query and subject lengths
# II) computing coverage and appending to BLAST output

# -- message on what this script does
cat <<EOF
-- this script computes query and subject coverage from BLASTP output
   and appends the coverage values as new columns to the BLAST output file
EOF

# -- default parameters
INPUT_FILE='output/blast_output/blast_results.tsv'
PROTEIN_INFO_FILE='data/protein_info_longest.csv' # -- only care for ID in col 1 and length in col 10 (also make sure to not use header)
COVERAGE_OUTPUT_FILE=$(basename "$INPUT_FILE" .tsv)_with_coverage.tsv

# -- arguments
while getopts "i:o:p:h" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) COVERAGE_OUTPUT_FILE="${OPTARG}" ;;
        p) PROTEIN_INFO_FILE="${OPTARG}" ;;
        h)
            echo "Usage: $0 [-i input_file] [-o coverage_output_filename] [-p protein_info_file]"
            echo "  -i    Input BLAST results file"
            echo "  -o    Output file name for coverage results"
            echo "  -p    Protein info file"
            echo "  -h    Show this help message"
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

# -- output dir will be the same as dir of input file
OUTPUT_DIR=$(dirname ${INPUT_FILE})
COVERAGE_OUTPUT_FILE=${OUTPUT_DIR}/${COVERAGE_OUTPUT_FILE}


echo "-- computing coverage for BLASTP results:"
echo "   INPUT : $INPUT_FILE"
echo "   PROTEIN INFO FILE : $PROTEIN_INFO_FILE"
echo "   OUTPUT: $COVERAGE_OUTPUT_FILE"

# ------------------------------------------------------------------

# I) preparing for covergae
# -- logic behind this code chunk:
#    1. create a file of format: protein_id length (without header) => $temp_protein_lengths
#    2. join BLAST output with $temp_protein_lengths on query id to get query lengths => $blast_with_query_lengths
#    3. join the result with $temp_protein_lengths on subject id to get subject lengths => $blast_with_subject_lengths
#    4. the final file has the format (in temp_blast_with_subject_lengths):
#       qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlength slength
# -- CREATING ALL TEMP FILES IN ONE FILE FOR DEBUGGING

mkdir -p "$OUTPUT_DIR/temp_files" 
temp_protein_lengths_file="${OUTPUT_DIR}/temp_files/temp_protein_lengths.txt"
blast_with_query_lengths_file="${OUTPUT_DIR}/temp_files/temp_blast_with_query_lengths.txt"
blast_with_subject_lengths_file="${OUTPUT_DIR}/temp_files/temp_blast_with_subject_lengths.txt"
header_file="${OUTPUT_DIR}/temp_files/temp_header.txt"
blast_with_coverage_file="${OUTPUT_DIR}/temp_files/temp_blast_with_coverage.txt"

# -- 1.
tail -n +2 "$PROTEIN_INFO_FILE" | awk -F, '{print $1, $10}' > "$temp_protein_lengths_file"
# -- 2.
join -1 1 -2 1 <(sort -k1,1 "$INPUT_FILE") <(sort -k1,1 "$temp_protein_lengths_file") | awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' > "$blast_with_query_lengths_file"
# -- 3.
join -1 2 -2 1 <(sort -k2,2 "$blast_with_query_lengths_file") <(sort -k1,1 "$temp_protein_lengths_file") | awk '{OFS="\t"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $14, $13}' > "$blast_with_subject_lengths_file"

echo "-- prepared BLAST output with query and subject lengths"

# ----------------------------------------------------------

# II) computing coverage
# -- the following steps:
#    1. compute query and subject coverage and append to BLAST output => temp_blast_with_coverage.txt
#    2. add header to final output file => $COVERAGE_OUTPUT_FILE (blast_results_with_coverage.tsv)

# -- function to compute coverage
compute_cov() {
    awk '{
        q_start = $7
        q_end   = $8
        s_start = $9
        s_end   = $10
        q_length = $13
        s_length = $14
        id = $3
        bit = $12

        q_cov = (q_end - q_start + 1) / q_length * 100
        s_cov = (s_end - s_start + 1) / s_length * 100

        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, q_cov, s_cov
    }' "$1"
}

# -- 1. compute coverage and add to BLAST output
compute_cov "${blast_with_subject_lengths_file}" > "$blast_with_coverage_file"
echo "-- computed coverage"
echo "   last 2 columns correspond to qcoverage and scoverage"

# -- 2. fix header based on qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcov scov
echo "qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore    qlength slength qcov	scov" > "$header_file"
cat "$header_file" "$blast_with_coverage_file" > "$COVERAGE_OUTPUT_FILE"

# ----------------------------------------------------------

# -- remove all temp files (keep for debugging)
rm -rf "${OUTPUT_DIR}/temp_files"
echo "-- removed temporary files"

echo "-- coverage computation completed, results saved to ${COVERAGE_OUTPUT_FILE}"

exit 0