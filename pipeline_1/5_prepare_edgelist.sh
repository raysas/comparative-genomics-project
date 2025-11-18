#!/bin/bash

# -- message on what this script does
cat <<EOF
-- this script generates an edgelist from the filtered BLASTP results
   the edgelist will be of the form: protein1 protein2 weight
   where weight here is bit score (consider making the option more general later)

EOF

# -- default parameters
INPUT_FILE='output/blast_filtered/filtered_blast_results_id50_qcov70_scov70_evalue1.tsv'
OUTPUT_DIR='output/similarity_edgelists'
WEIGHT_COLUMN_INDEX=12  # -- bit score column in BLAST output

# -- get arguments
while getopts "i:o:h" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        w) WEIGHT_COLUMN_INDEX="${OPTARG}" ;;
        h)
            echo "Usage: $0 [-i input_file] [-o OUTPUT_DIR] [-w weight_column_index]"
            echo "  -i    Input filtered BLAST results file"
            echo "  -o    Output edgelist file"
            echo "  -w    Column index for weight in BLAST output"
            echo "  -h    Show this help message"
            exit 0
            ;;
        *)
            echo "Invalid option. Use -h for help."
            exit 1
            ;;
    esac
done

# -- if weight index in input file is not a number, say sorry the weight should be a number and set to bit score column by default

LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR" 
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
echo "Command: $0 $*"

exec > >(tee -i "$LOG_FILE") 2>&1

OUTPUT_FILE=$OUTPUT_DIR/$(basename "${INPUT_FILE%.tsv}_wcol${WEIGHT_COLUMN_INDEX}.txt")

echo "-- generating edgelist from filtered BLASTP results:"
echo "   INPUT : $INPUT_FILE"
echo "   OUTPUT: $OUTPUT_FILE"
echo "   WEIGHT COLUMN INDEX: $WEIGHT_COLUMN_INDEX"

tail -n +2 "$INPUT_FILE" |awk -v weight_col="$WEIGHT_COLUMN_INDEX" 'BEGIN {
    if (weight_col !~ /^[0-9]+$/) {
        print "Error: Weight column index must be a number. Setting to default bit score column (12)"
        weight_col = 12
    }
}
{
    print $1, $2, $weight_col
}' > "$OUTPUT_FILE"

echo "-- edgelist written to $OUTPUT_FILE"

exit 0