#!/usr/bin/env bash
set -euo pipefail

ctl_dir="../output/ks/ctl"
results_dir="../output/ks/results"

start=1
end=7221

for ((fam_num=start; fam_num<=end; fam_num++)); do
    fam_base="family${fam_num}"
    fam_ctl_dir="$ctl_dir/$fam_base"
    fam_results_dir="$results_dir/$fam_base"

    echo "▶ Processing $fam_base"
    mkdir -p "$fam_results_dir"

    shopt -s nullglob
    for ctl in "$fam_ctl_dir"/*.yn00.ctl; do
        pair_base=$(basename "$ctl" .yn00.ctl)
        pair_dir="$fam_results_dir/$pair_base"
        mkdir -p "$pair_dir"

        echo "   ▶ Running yn00 for $pair_base"
        yn00 "$ctl"

      
        mv yn rst rst1 rub 2YN.t 2YN.dN 2YN.dS "$pair_dir/" 2>/dev/null || true
    done
    shopt -u nullglob
done

echo "Ka/Ks results saved in subfolders for each gene pair under $results_dir"
