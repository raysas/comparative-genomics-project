#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    clusters proteins into families based on filtered BLASTP results
#    using connected components from a graph representation of the BLAST results
#    will be using the MCL algorithm for clustering
# !! need to fix merged ids: KRH29797KRH39445
# !! need to format to tsv: geneName\tfamily
#
#    so 2 formats of output:
#    1) .txt -mcl output format: each line is a family with space-separated gene
#    2) .tsv - 2 columns: geneName and familyID
#
# -- Usage:
#    bash ./pipeline_1/6_cluster_families.sh [-i INPUT_FILE] [-o OUTPUT_DIR] [-h]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/6_cluster_families.sh -i "output/similarity_edgelists/filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv" -o "output/clusters"
# --------------------------------------------------------------------

# -- message on what this script does
cat <<EOF
-- this script clusters proteins into families based on filtered BLASTP results
    using connected components from a graph representation of the BLAST results
    will be using the MCL algorithm for clustering
EOF

# -- default parameters
INPUT_FILE='output/similarity_edgelists/filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv'
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
TSV_OUTPUT_FILE="${OUTPUT_FILE%.txt}.tsv"

# -- this is the default command for mcl ran on galaxy
mcl "$INPUT_FILE" -I "$MCL_INFLATION" --abc -V all --discard-loops="$DISCARD_LOOPS" -c 1.0 -P "$PRUNING_THRESHOLD" -S "$SELECT_DOWN_TO" -R "$RECOVER" -pct "$PCT" -o "$OUTPUT_FILE"
echo "-- clustering completed"

# -- 2 things to do:
# * fix mixed ids in output (merged ids) e.g. KRH29797KRH39445 should be KRH29797\tKRH39445
# * map all ids in the same line to the same familyID

echo "-- converting MCL output into TSV: $TSV_OUTPUT_FILE"
echo "-- fixing merged ids and mapping to family IDs"

python3 <<EOF
import re

input_file = "$OUTPUT_FILE"
tsv_file = "$TSV_OUTPUT_FILE"

# --- step 1: fix merged IDs and rewrite the MCL file in place ---
cleaned_lines = []
with open(input_file) as f:
    for line in f:
        line = line.rstrip()
        # Fix merged IDs: insert space between digit → letter, e.g., "123A"
        fixed = re.sub(r"(?<=\d)(?=[A-Za-z])", " ", line)
        cleaned_lines.append(fixed)

# overwriting the MCL output file with cleaned lines
with open(input_file, "w") as f:
    for line in cleaned_lines:
        f.write(line + "\n")

# --- step 2: generate TSV from the now-corrected file ---
with open(input_file) as f, open(tsv_file, "w") as out:
    out.write("geneName\tfamily\n")
    family_id = 1
    for line in f:
        genes = [g for g in line.split() if g]
        for g in genes:
            out.write(f"{g}\t{family_id}\n")
        family_id += 1
EOF

echo "-- TSV conversion complete: $TSV_OUTPUT_FILE"

# -- transform to tsv format: geneName\tfamilyID

exit 0