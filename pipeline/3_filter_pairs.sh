#!/bin/sh

# -- meaasge what this script does

cat <<EOF
-- this script filters BLASTP results based on identity/coverge/score thresholds
   and outputs the filtered results to a new file
   can specify the features to remove from arguments, through the index of the column in the BLAST output with coverage added
EOF

# -- default parameters
INPUT_FILE='output/blast_output/blast_results_with_coverage.tsv'
OUTPUT_FILE='output/blast_output/filtered_blast_results.tsv'

ID_THRESHOLD=30
Q_COV_THRESHOLD=0.5
S_COV_THRESHOLD=0.5
BIT_SCORE_THRESHOLD=NULL  # -- no filtering on bit score by default
E_VALUE_THRESHOLD=NULL  # -- no filtering on e-value by default

LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1

# -- get input args from command line
while getopts i:o:id:qcov:scov:bit:evalue:h flag
do
    case "${flag}" in
        i) INPUT_FILE=${OPTARG};;
        o) OUTPUT_FILE=${OPTARG};;
        id) ID_THRESHOLD=${OPTARG};;
        qcov) Q_COV_THRESHOLD=${OPTARG};;
        scov) S_COV_THRESHOLD=${OPTARG};;
        bit) BIT_SCORE_THRESHOLD=${OPTARG};;
        evalue) E_VALUE_THRESHOLD=${OPTARG};;
        h) echo "Usage: $0 [-i input_file] [-o output_file] [-id identity_threshold] [-qcov query_coverage_threshold] [-scov subject_coverage_threshold] [-bit bit_score_threshold] [-evalue e_value_threshold]"
           exit 0;;
    esac
done

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

# -- helper function to remove lines that are less than or greater than
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

# -- for identity its remove less than threshold
filter_by_feature "$INPUT_FILE" 3 "$ID_THRESHOLD" "less" "$OUTPUT_FILE.tmp1"
# -- for query coverage
filter_by_feature "$OUTPUT_FILE.tmp1" 13 "$Q_COV_THRESHOLD" "less" "$OUTPUT_FILE.tmp2"
# -- for subject coverage
filter_by_feature "$OUTPUT_FILE.tmp2" 14 "$S_COV_THRESHOLD" "less" "$OUTPUT_FILE.tmp3"
# -- for bit score
filter_by_feature "$OUTPUT_FILE.tmp3" 12 "$BIT_SCORE_THRESHOLD" "less" "$OUTPUT_FILE.tmp4"
# -- for e-value
filter_by_feature "$OUTPUT_FILE.tmp4" 11 "$E_VALUE_THRESHOLD" "greater" "$OUTPUT_FILE.tmp5"

mv "$OUTPUT_FILE.tmp5" "$OUTPUT_FILE"
rm "$OUTPUT_FILE.tmp1" "$OUTPUT_FILE.tmp2" "$OUTPUT_FILE.tmp3" "$OUTPUT_FILE.tmp4" "$OUTPUT_FILE.tmp5"

echo "-- filtering completed! filtered results saved to $OUTPUT_FILE"

# -- end of script
echo "-- end of script $0"
exit 0