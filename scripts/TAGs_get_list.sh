#!/bin/bash

input_dir="output/duplication_classes/TAGs/"
output_dir="output/gene_lists/TAGs/spacer_based/"

if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

echo "-- input directory: $input_dir"
# iterate over subdirectories safely (nullglob avoids literal pattern if no matches)
shopt -s nullglob
for subdir in "${input_dir}"*/; do
    # ensure it's a directory
    [ -d "$subdir" ] || continue
    input_file="${subdir}TAGs_1.tsv"
    output_file="${output_dir}TAGs_$(basename "$subdir").txt"
    echo "-- Processing $input_file"
    if [ -f "$input_file" ]; then
        # assume CSV with commas; use -F, to split on comma
        awk 'NR>1 && $6>0 {print $4}' "$input_file" > "$output_file"
        echo "-- IDs saved to $output_file"
    else
        echo "-- Warning: input file not found: $input_file (skipping)"
    fi
done
shopt -u nullglob
echo "-- can find them in $output_dir"
