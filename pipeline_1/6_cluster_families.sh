#!/bin/bash

# -- message on what this script does
cat <<EOF
-- this script clusters proteins into families based on filtered BLASTP results
    using connected components from a graph representation of the BLAST results
    will be using the MCL algorithm for clustering
EOF

# -- default parameters
INPUT_FILE='output/similarity_edgelists/cluster_input_id50_qcov70_scov70_evalue1_wcol12.tsv'
OUTPUT_DIR='output/clusters'
MCL_INFLATION=2.0
DISCARD_LOOPS='y'
PRUNING_THRESHOLD=4000
SELECT_DOWN_TO=500
RECOVER=600
PCT=90
C_FOLD=1.0

# -- get arguments
while getopts "i:o:h" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        I) MCL_INFLATION="${OPTARG}" ;;
        P) PRUNING_THRESHOLD="${OPTARG}" ;;
        S) SELECT_DOWN_TO="${OPTARG}" ;;
        R) RECOVER="${OPTARG}" ;;
        discard_loops) DISCARD_LOOPS="${OPTARG}" ;;
        pct) PCT="${OPTARG}" ;;
        c) C_FOLD="${OPTARG}" ;;
        h)
            echo "Usage: $0 [-i input_file] [-o output_dir]"
            echo "  -i    Input edgelist file for clustering"
            echo "  -o    Output file for protein families"
            echo "  -h    Show this help message"
            echo "  -c <num>           increase loop-weights <num>-fold (default: $C_FOLD)"
            echo "  -discard_loops y|n   discard self-loops (default: $DISCARD_LOOPS)"
            echo "  -I <num>          MCL inflation parameter (default: $MCL_INFLATION)"
            echo "  -P <num>          MCL (inverted) rigid pruning threshold (cf -z) (default: $PRUNING_THRESHOLD)"
            echo "  -S <num>          MCL select down to <int> entries if needed"
            echo "  -R <num>          MCL recover to maximally <int> entries if needed"
            echo "  -pct <num>        MCL try recovery if mass is less than <pct>"
            echo "for more info see mcl --help"
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

echo "-- clustering proteins into families using MCL with parameters:"
echo "   INPUT : $INPUT_FILE"
echo "   OUTPUT: $OUTPUT_DIR"
echo "   MCL Inflation: $MCL_INFLATION"
echo "   Discard Loops: $DISCARD_LOOPS"
echo "   Pruning Threshold: $PRUNING_THRESHOLD"
echo "   Select Down To: $SELECT_DOWN_TO"
echo "   Recover: $RECOVER"
echo "   PCT: $PCT" 
echo "   C-Fold: $C_FOLD"

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    echo "-- created output directory: $OUTPUT_DIR"
fi

# -- rename file based on input file
OUTPUT_FILE="${OUTPUT_DIR}/protein_families_$(basename "${INPUT_FILE%.tsv}").txt"

# -- this is the default command for mcl ran on galaxy
mcl "$INPUT_FILE" -I "$MCL_INFLATION" --abc -V all --discard-loops="$DISCARD_LOOPS" -c 1.0 -P "$PRUNING_THRESHOLD" -S "$SELECT_DOWN_TO" -R "$RECOVER" -pct "$PCT" -o "$OUTPUT_FILE"

echo "-- clustering completed: Output written to $OUTPUT_FILE"

exit 0