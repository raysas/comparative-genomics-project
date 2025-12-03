#!/usr/bin/env bash
set -euo pipefail

tsv="../output/clusters/protein_families_network_50families_max5.tsv"
prot_db="../data/peptides.fa"
cds_db="../data/cds.fa"
out_dir="../output/pairwise"

mkdir -p "$out_dir"

# Lặp qua từng family
awk 'NR>1{print $2}' "$tsv" | sort | uniq | while read fam; do
    fam_clean=$(echo "$fam" | tr -cd '0-9')   # chỉ giữ số
    echo "▶ Processing family $fam_clean"

    ids=$(awk -v f="$fam" 'NR>1 && $2==f {print $1}' "$tsv")

    for id1 in $ids; do
        for id2 in $ids; do
            if [[ "$id1" < "$id2" ]]; then
                pair="family${fam_clean}_${id1}_${id2}"
                echo "  ▶ Pair $id1 vs $id2"

                # Extract protein sequences (match header chứa ID)
                awk -v id1="$id1" -v id2="$id2" \
                    '/^>/ {keep=($0 ~ id1 || $0 ~ id2)} keep' "$prot_db" \
                    > "$out_dir/${pair}.prot.fst"

                # Extract CDS sequences
                awk -v id1="$id1" -v id2="$id2" \
                    '/^>/ {keep=($0 ~ id1 || $0 ~ id2)} keep' "$cds_db" \
                    > "$out_dir/${pair}.cds.fst"

                # Skip if missing sequences
                if [[ ! -s "$out_dir/${pair}.prot.fst" || ! -s "$out_dir/${pair}.cds.fst" ]]; then
                    echo "⚠️ Skipping $pair (missing sequences)"
                    continue
                fi

                # Step 1: protein alignment
                clustalw2 -quiet -align -infile="$out_dir/${pair}.prot.fst" \
                          -outfile="$out_dir/${pair}.prot.aln"

                # Step 2: CDS alignment
                ./pal2nal.pl "$out_dir/${pair}.prot.aln" "$out_dir/${pair}.cds.fst" \
                             -output paml > "$out_dir/${pair}.cds.phy"

                # Step 3: control file
                awk -v file="$out_dir/${pair}.cds.phy" \
                    '{gsub("XXXXX",file); print $0}' yn00.ctl_master \
                    > "$out_dir/${pair}.yn00.ctl"

                # Step 4: run yn00
                yn00 "$out_dir/${pair}.yn00.ctl"
                mv yn "$out_dir/${pair}.yn"
            fi
        done
    done
done

echo "✅ Pairwise Ka/Ks results saved in $out_dir"
