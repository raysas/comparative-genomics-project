#!/usr/bin/env bash
set -euo pipefail

# This script runs yn00 for Ka/Ks calculation.
# Currently hardcoded to use the "test" dataset.
# To switch to "full" dataset, change the paths below.

ctl_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/test/ctl"
results_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/test/results"

# Select the range of families to run
start=1     # starting family
end=20      # ending family

mkdir -p "$results_dir"

# Loop from start to end
for ((fam_num=start; fam_num<=end; fam_num++)); do
    fam_base="family${fam_num}"
    ctl="$ctl_dir/${fam_base}.yn00.ctl"

    if [[ -f "$ctl" ]]; then
        echo "▶ Running yn00 for $fam_base"
        yn00 "$ctl"
    else
        echo "Skipping $fam_base (no control file found)"
    fi
done

echo "Ka/Ks results saved in $results_dir for families $start to $end"
