#!/bin/bash

input_file="output/duplication_classes/proximal/proximal_500kb.csv"
output_file="output/gene_lists/proximal.txt"
tag_list="output/gene_lists/TAGs.txt"

# Skip header and print peptide_id if proximal > 0
awk 'NR>1 && $6>0 {print $4}' "$input_file" > "$output_file"
echo "-- extracted proximal gene IDs from $input_file"

# -- remove intersection with TAGs
grep -v -F -f "$tag_list" "$output_file" > "${output_file}.tmp"
mv "${output_file}.tmp" "$output_file"
echo "-- removed proximal gene IDs that are also TAGs from $tag_list"

echo "-- IDs saved to $output_file"
