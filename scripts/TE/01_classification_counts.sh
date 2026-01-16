#!/bin/bash
# Script Name: 01_classification_counts.sh
# Purpose: Generate counts of unique TE Classifications from the GFF file.

# Define Input and Output directories to keep things organized
INPUT_FILE="../data/TEAnnotationFinal.gff3"
OUTPUT_FILE="../output/TE_classification_counts.txt"

# Ensure output directory exists
mkdir -p ../output

echo "Analyzing TE classifications in $INPUT_FILE..."

# ---------------------------------------------------------
# The Analysis Pipeline
# ---------------------------------------------------------
# 1. grep -v "#": Reads the file and removes header lines (comments).
# 2. cut -f3    : Extracts Column 3 (The Feature/Classification type).
# 3. sort       : Sorts names alphabetically (required for uniq to work).
# 4. uniq -c    : Collapses identical lines and counts them.
# 5. sort -nr   : Sorts the final counts numerically (n) in reverse (r).
# ---------------------------------------------------------

grep -v "#" "$INPUT_FILE" | cut -f3 | sort | uniq -c | sort -nr > "$OUTPUT_FILE"

echo "Analysis complete. Displaying results:"
echo "------------------------------------------------"

# Display the results to the terminal
cat "$OUTPUT_FILE"