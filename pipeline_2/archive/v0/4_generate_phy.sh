#!/usr/bin/env bash
set -euo pipefail

cds_dir="../output/ks/cds"
prot_align_dir="../output/ks/alignments"
phylip_dir="../output/ks/phylip"

# Select the range of families to run
start=4666
end=4667

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
        # Only process files ending with _prot.aln
        [[ "$aln" == *_prot.aln ]] || continue

        pair_base=$(basename "$aln" _prot.aln)
        cds_file="$cds_dir/$fam_base/${pair_base}.fa"
        out_file="$fam_out/${pair_base}.phy"

        echo "   ▶ Running pal2nal for $pair_base"

        # Remove CLUSTAL header line
        clean_aln=$(mktemp)
        sed '/^CLUSTAL/d' "$aln" > "$clean_aln"

        # Run pal2nal
        ./pal2nal.pl "$clean_aln" "$cds_file" -output paml > "$out_file"
        echo " -- Generated PHYLIP file: $out_file"

        rm "$clean_aln"
    done
done

echo "Pairwise CDS alignments saved in $phylip_dir/family$start to family$end/"
