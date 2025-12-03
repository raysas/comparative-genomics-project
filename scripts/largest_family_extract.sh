#!/bin/bash

# -- script aimed to extract the largest protein family from a given dataset
# -- the input dataset here is a file in output/clusters/
# -- usually files of the form:
#   protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv
#   protein_families_filtered_blast_results_id30_qcov50_scov50_evalue1e-5_wcol12_network.tsv


# -- 2 options:
# -- either provide -id -cov -evalue(optional) and will get the file myself
# -- or let user input path of the file

#!/bin/bash

# ------------------------------------------------------------
# Script to extract the largest protein family from a dataset
# Input file can be:
#   1) Automatically determined using -id -cov [-evalue]
#   2) Directly provided using -f /path/to/file
# ------------------------------------------------------------

# default values
ID=""
COV=""
EVALUE=""
FILE=""
OUTPUT_DIR="output/gene_lists/largest_family/"

usage() {
    echo "Usage:"
    echo "  Option 1 (auto-detect file):"
    echo "    $0 -id <identity> -cov <coverage> [-evalue <evalue>]"
    echo
    echo "  Option 2 (manual file path):"
    echo "    $0 -f <path/to/file>"
    exit 1
}

# ------------------------
# Parse command-line args
# ------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -id)
            ID="$2"
            shift 2
            ;;
        -cov)
            COV="$2"
            shift 2
            ;;
        -evalue)
            EVALUE="$2"
            shift 2
            ;;
        -f)
            FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done



if [[ -n "$FILE" ]]; then
    # user provided a file path, use it
    infile="$FILE"
    outfile="${OUTPUT_DIR}largest_family_$(basename "$FILE" | sed 's/protein_families_filtered_blast_results//; s/_network.tsv//').txt"
else
    # user should provide id and cov to auto-detect file
    if [[ -z "$ID" || -z "$COV" ]]; then
        echo "Error: You must provide either -f <file> or -id and -cov (and optionally -evalue)."
        usage
    fi
    if [[ -n "$EVALUE" ]]; then
        infile="output/clusters/protein_families_filtered_blast_results_id${ID}_qcov${COV}_scov${COV}_evalue${EVALUE}_wcol12_network.tsv"
        outfile="${OUTPUT_DIR}largest_family_id${ID}_cov${COV}_evalue${EVALUE}.txt"
    else
        infile="output/clusters/protein_families_filtered_blast_results_id${ID}_qcov${COV}_scov${COV}_wcol12_network.tsv"
        outfile="${OUTPUT_DIR}largest_family_id${ID}_cov${COV}.txt"
    fi
fi

if [ ! -f "$infile" ]; then
    echo "Error: input file not found: $infile" >&2
    exit 1
fi

read count family_id <<<$(awk '{count[$2]++} END{for (f in count) print count[f], f}' "$infile" | sort -nr | head -1)
echo "-- largest family ID: $family_id with $count members"

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi
awk -v fam_id="$family_id" '$2 == fam_id {print $1}' "$infile" > "$outfile"
echo "-- largest family members written to: $outfile"