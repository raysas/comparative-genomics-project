#!/usr/bin/env bash
set -euo pipefail

phy_dir="../output/ks/phylip"
ctl_dir="../output/ks/ctl"

start=4666
end=4667

mkdir -p "$ctl_dir"

for ((fam_num=start; fam_num<=end; fam_num++)); do
    fam_base="family${fam_num}"
    fam_phy_dir="$phy_dir/$fam_base"
    fam_ctl_dir="$ctl_dir/$fam_base"
    mkdir -p "$fam_ctl_dir"

    echo "▶ Processing $fam_base"

    shopt -s nullglob
    for phy in "$fam_phy_dir"/*.phy; do
        pair_base=$(basename "$phy" .phy)
        ctl_file="$fam_ctl_dir/${pair_base}.yn00.ctl"

        echo "   ▶ Generating control file for $pair_base"
        awk -v file="$phy" '{gsub("XXXXX",file); print $0}' yn00.ctl_master > "$ctl_file"
    done
    shopt -u nullglob
done

echo "Control files saved in $ctl_dir for families $start to $end"
