#!/usr/bin/env bash
set -euo pipefail

# This script runs pal2nal.pl to generate CDS alignments in PHYLIP format.
# It is configured for the "test" dataset. Change paths for "full" dataset.

cds_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/full/families"
prot_align_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/full/alignments"
phylip_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/full/phylip"

# Select the range of families to run
start=4666     # starting family
end=4667    # ending family

# Ensure the output folder exists
mkdir -p "$phylip_dir"

# Loop from start to end
for ((fam_num=start; fam_num<=end; fam_num++)); do
    fam_base="family${fam_num}"
    fam_dir="$prot_align_dir/$fam_base"
    fam_out="$phylip_dir/$fam_base"
    mkdir -p "$fam_out"

    echo "▶ Processing $fam_dir -> $fam_out"

    # Loop through all pairwise protein alignment files in this family
    for aln in "$fam_dir"/*.aln; do
        pair_base=$(basename "$aln" .aln)
        cds_file="$cds_dir/$fam_base/${pair_base}.fa"
        out_file="$fam_out/${pair_base}.phy"

        echo "   ▶ Running pal2nal for $pair_base"

        # Remove CLUSTAL header line from alignment file
        clean_aln=$(mktemp)
        sed '/^CLUSTAL/d' "$aln" > "$clean_aln"

        # Run pal2nal with the cleaned alignment and corresponding CDS file
        ./pal2nal.pl "$clean_aln" "$cds_file" -output paml > "$out_file"

        # Delete temporary cleaned alignment file
        rm "$clean_aln"
    done
done

echo "Pairwise CDS alignments saved in $phylip_dir/family$start to family$end/"
