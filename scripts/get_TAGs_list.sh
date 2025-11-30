#!/bin/bash

input_file="output/duplication_classes/TAGs/TAGs_1.tsv"
output_file="output/duplication_classes/TAGs.txt"

# Skip header and print peptide_id if TAG > 0
awk 'NR>1 && $6>0 {print $4}' "$input_file" > "$output_file"

echo "-- IDs saved to $output_file"
