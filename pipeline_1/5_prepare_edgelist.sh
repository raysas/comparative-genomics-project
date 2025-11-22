#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    generates an edgelist from the filtered BLASTP results
#    the edgelist will be of the form: protein1 protein2 weight
#    where weight here is bit score (consider making the option more general later)
#
#  !SOME ISSUES WITH THE SIMILARITY NETWORK:
#    we have symmetric edges (p1,p2,w) and (p2,p1,w) that needs to be shortned to one edge
#    we always have a self loop (best match to p1 is p1): needs to be removed
#    we have multiple edges between some pairs (p1,p2) with different weights (due to different alignments): need to keep only the highest weight
#    * will create temp files during filtreation for testing and troubleshooting
# -- Usage:
#    bash ./pipeline_1/5_prepare_edgelist.sh [-i INPUT_FILE] [-o OUTPUT_DIR] [-w WEIGHT_COLUMN_INDEX] [-h]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/5_prepare_edgelist.sh -i "output/blast_filtered/filtered_blast_results_id30_qcov50_scov50.tsv" -o "output/similarity_edgelists" -w 12
# --------------------------------------------------------------------
###########################################################################

# ---------------------------------------------------------------------
# prepare variables, get arguments, set up logging
# ---------------------------------------------------------------------

# -- message on what this script does
cat <<EOF
-- this script generates an edgelist from the filtered BLASTP results
   the edgelist will be of the form: protein1 protein2 weight
   where weight here is bit score (consider making the option more general later)

EOF

# -- default parameters
INPUT_FILE='output/blast_filtered/filtered_blast_results_id30_qcov50_scov50.tsv'
OUTPUT_DIR='output/similarity_edgelists'
WEIGHT_COLUMN_INDEX=12  # -- bit score column in BLAST output

# -- get arguments
while getopts "i:o:w:h" flag; do
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

# -- check input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: input file $INPUT_FILE not found" >&2
    exit 1
fi
# -- ensure output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    echo "-- created output directory $OUTPUT_DIR"
fi


OUTPUT_FILE=$OUTPUT_DIR/$(basename "${INPUT_FILE%.tsv}_wcol${WEIGHT_COLUMN_INDEX}_network.tsv")

echo " -- generating edgelist from filtered BLASTP results:"
echo "   INPUT : $INPUT_FILE"
echo "   OUTPUT: $OUTPUT_FILE"
echo "   WEIGHT COLUMN INDEX: $WEIGHT_COLUMN_INDEX"
echo " -- starting with $(wc -l < "$INPUT_FILE") lines in input file"

# -------------------------------------------------------------------------
# -- main script logic
# -------------------------------------------------------------------------

# -- re-order columns to have first id < second id (lexico-alphabetical order)

# -- extract relevant columns: protein1, protein2, weight

tail -n +2 "$INPUT_FILE" | awk -v weight_col="$WEIGHT_COLUMN_INDEX" 'BEGIN {
    if (weight_col !~ /^[0-9]+$/) {
        print "Error: Weight column index must be a number. Setting to default bit score column (12)" > "/dev/stderr"
        weight_col = 12
    }
}
{
    print $1, $2, $(weight_col)
}' > "$OUTPUT_FILE"

# -- reorder to have first id < second id (lexico-alphabetical order)
# -- !! made sure to have the real query and subject ids in a row before the swapped ones (by sorting first)
sort -d ${OUTPUT_FILE} | awk '
{
    if ($1 <= $2) {
        print $0
    } else {
        # --swap $1 and $2, keep rest untouched
        temp = $1
        $1 = $2
        $2 = temp
        print $0
    }
}
' > "${OUTPUT_FILE}.tmp"
echo "-- reordered columns to ensure first id <= second id lexicographically"

# -- 1.remove self loops (when col1=col2)
awk '$1 != $2' "$OUTPUT_FILE" > "${OUTPUT_FILE}.tmp"
mv "${OUTPUT_FILE}.tmp" "$OUTPUT_FILE"
rm -f "${OUTPUT_FILE}.tmp"
echo "-- self-loops removed"
echo "-- $(wc -l < "$OUTPUT_FILE") lines remain after removing self-loops"

# -- 2.some duplicates of col1+col2, will keep only the highest weight (ALREADY SORTED)
cat "$OUTPUT_FILE" | awk '{
    key = $1 "_" $2
    if (!(key in seen)) {
        print $0
        seen[key] = 1
    }
}' > "${OUTPUT_FILE}.tmp" 

mv "${OUTPUT_FILE}.tmp" "$OUTPUT_FILE"
rm -f "${OUTPUT_FILE}.tmp"
echo "-- duplicate edges removed, keeping highest weight"
echo "-- $(wc -l < "$OUTPUT_FILE") lines remain after removing duplicate edges with different weights (multiple hits - keeping best hit)"

echo "-- edgelist written to $OUTPUT_FILE"

exit 0