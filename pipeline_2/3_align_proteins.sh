#!/usr/bin/env bash
set -euo pipefail

# It is configured for the "test" dataset. Change paths for "full" dataset.
input_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/full/families"
output_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/full/alignments"

# Select the range of families to run
start=4666    # starting family
end=4667      # ending family

# Ensure the output folder exists
mkdir -p "$output_dir"

# Loop from start to end
for ((fam_num=start; fam_num<=end; fam_num++)); do
    fam_dir="$input_dir/family${fam_num}"
    fam_out="$output_dir/family${fam_num}"
    mkdir -p "$fam_out"

    echo "▶ Processing $fam_dir -> $fam_out"

    # Loop through all pairwise FASTA files in this family
    for pair in "$fam_dir"/*.fa; do
        pair_base=$(basename "$pair" .fa)
        out_file="$fam_out/${pair_base}.aln"
        echo "   ▶ Aligning $pair -> $out_file"
        clustalw2 -quiet -align -infile="$pair" -outfile="$out_file"
    done
done

echo "Pairwise alignments done for families $start to $end, results in $output_dir/family*/"