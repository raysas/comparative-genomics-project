#!/bin/bash

# -- meaasge what this script does

cat <<EOF
-- this script filters BLASTP results based on identity/coverge/score thresholds
   and outputs the filtered results to a new file
   can specify the features to remove from arguments, through the index of the column in the BLAST output with coverage added
EOF

# -- default parameters
INPUT_FILE='output/blast_output/blast_results_with_coverage.tsv'
OUTPUT_FILE='output/blast_output/filtered_blast_results'

ID_THRESHOLD=30
Q_COV_THRESHOLD=50
S_COV_THRESHOLD=50
BIT_SCORE_THRESHOLD=NULL  # -- no filtering on bit score by default
E_VALUE_THRESHOLD=NULL  # -- no filtering on e-value by default

# -- get input args from command line (changed a bit to fix an error)
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) INPUT_FILE="$2"; shift 2;;
        -o) OUTPUT_FILE="$2"; shift 2;;
        -id) ID_THRESHOLD="$2"; shift 2;;
        -qcov) Q_COV_THRESHOLD="$2"; shift 2;;
        -scov) S_COV_THRESHOLD="$2"; shift 2;;
        -bit) BIT_SCORE_THRESHOLD="$2"; shift 2;;
        -evalue) E_VALUE_THRESHOLD="$2"; shift 2;;
        -h|--help)
            echo "Usage: $0 [-i input_file] [-o output_file] [-id identity_threshold] [-qcov query_coverage_threshold] [-scov subject_coverage_threshold] [-bit bit_score_threshold] [-evalue e_value_threshold]"
            echo "  -i        Input BLAST results file with coverage"
            echo "  -o        Output file name for filtered results"
            echo "  -id       Identity threshold (default: $ID_THRESHOLD)"
            echo "  -qcov     Query coverage threshold (default: $Q_COV_THRESHOLD)"
            echo "  -scov     Subject coverage threshold (default: $S_COV_THRESHOLD)"
            echo "  -bit      Bit score threshold (default: no filtering)"
            echo "  -evalue   E-value threshold (default: no filtering)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"; exit 1;;
    esac
done

LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"


echo " -- this file parameters:"
echo "   INPUT FILE : $INPUT_FILE"
echo "   OUTPUT FILE: $OUTPUT_FILE"
echo " -- filtering BLASTP results in $INPUT_FILE with thresholds:"
echo "   Identity >= $ID_THRESHOLD"
echo "   Query Coverage >= $Q_COV_THRESHOLD"
echo "   Subject Coverage >= $S_COV_THRESHOLD"
if [ "$BIT_SCORE_THRESHOLD" != "NULL" ]; then
    echo "   Bit Score >= $BIT_SCORE_THRESHOLD"
else
    echo "   No filtering on Bit Score"
fi
if [ "$E_VALUE_THRESHOLD" != "NULL" ]; then
    echo "   E-value <= $E_VALUE_THRESHOLD"
else
    echo "   No filtering on E-value"
fi

# -- put thresholds in output filename
OUTPUT_FILE="${OUTPUT_FILE}_id${ID_THRESHOLD}_qcov${Q_COV_THRESHOLD}_scov${S_COV_THRESHOLD}"
if [ "$BIT_SCORE_THRESHOLD" != "NULL" ]; then
    OUTPUT_FILE="${OUTPUT_FILE}_bit${BIT_SCORE_THRESHOLD}"
fi
if [ "$E_VALUE_THRESHOLD" != "NULL" ]; then
    OUTPUT_FILE="${OUTPUT_FILE}_evalue${E_VALUE_THRESHOLD}"
fi
OUTPUT_FILE="${OUTPUT_FILE}.tsv"

################################################################################
# ------------------------------ old code to delete ---------------------------
# -- function that takes column index and removes lines with value below threshold
# remove_less_than() {
#     input_file="$1"
#     column_index="$2"
#     threshold="$3"
#     output_file="$4"

#     awk -v col_idx="$column_index" -v thresh="$threshold" '{
#         if ($col_idx >= thresh) {
#             print $0
#         }
#     }' "$input_file" > "$output_file"
# }
# remove_greater_than() {
#     input_file="$1"
#     column_index="$2"
#     threshold="$3"
#     output_file="$4"

#     awk -v col_idx="$column_index" -v thresh="$threshold" '{
#         if ($col_idx <= thresh) {
#             print $0
#         }
#     }' "$input_file" > "$output_file"
# }
# ------------------------------ end old code to delete -----------------------
################################################################################

# -- helper function to remove lines that are less than or greater than
# -- how to use it:
#       filter_by_feature input_file column_index threshold comparison[less/greater] output_file
filter_by_feature() {
    input_file="$1"
    column_index="$2"
    threshold="$3"
    comparison="$4"  # -- "less" or "greater"
    output_file="$5"

    # -- if threshold is NULL, skip filtering
    if [ "$threshold" = "NULL" ]; then
        cp "$input_file" "$output_file"
        echo "-- skipping filtering for column $column_index as threshold is NULL"
        return
    fi

    echo "-- filtering column $column_index with threshold $threshold (remoing values $comparison than threshold)"

    if [ "$comparison" = "less" ]; then
        awk -v col_idx="$column_index" -v thresh="$threshold" '{
            if ($col_idx >= thresh) {
                print $0
            }
        }' "$input_file" > "$output_file"
    elif [ "$comparison" = "greater" ]; then
        awk -v col_idx="$column_index" -v thresh="$threshold" '{
            if ($col_idx <= thresh) {
                print $0
            }
        }' "$input_file" > "$output_file"
    else
        echo "Error: comparison must be 'less' or 'greater'"
        exit 1
    fi
}

# -- temp folder
TEMP_DIR="$(dirname "$OUTPUT_FILE")/temp/"
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p "$TEMP_DIR"
fi

# -- for identity its remove less than threshold
filter_by_feature "$INPUT_FILE" 3 "$ID_THRESHOLD" "less" $TEMP_DIR/by_id.tsv
# -- for query coverage
filter_by_feature  $TEMP_DIR/by_id.tsv 13 "$Q_COV_THRESHOLD" "less" "$TEMP_DIR/by_qcov.tsv"
# -- for subject coverage
filter_by_feature "$TEMP_DIR/by_qcov.tsv" 14 "$S_COV_THRESHOLD" "less" "$TEMP_DIR/by_scov.tsv"
# -- for bit score
filter_by_feature "$TEMP_DIR/by_scov.tsv" 12 "$BIT_SCORE_THRESHOLD" "less" "$TEMP_DIR/by_bitscore.tsv"
# -- for e-value
filter_by_feature "$TEMP_DIR/by_bitscore.tsv" 11 "$E_VALUE_THRESHOLD" "greater" "$TEMP_DIR/by_evalue.tsv"


echo "-- filtering completed! filtered results saved to $OUTPUT_FILE"
cp "$TEMP_DIR/by_evalue.tsv" "$OUTPUT_FILE"

# -- remove temp folder
rm -rf "$TEMP_DIR"
echo "-- removed temporary files in $TEMP_DIR"

exit 0