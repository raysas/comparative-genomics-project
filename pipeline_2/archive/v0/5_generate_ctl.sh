#!/usr/bin/env bash
set -euo pipefail

# This script runs pal2nal.pl to generate CDS alignments in PHYLIP format.
# Currently configured to use the "test" dataset. Change paths for "full" dataset.

phy_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/test/phylip"
ctl_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks/test/ctl"

# Select the range of families to run
start=1     # starting family
end=20      # ending family

# Ensure output folder exists
mkdir -p "$ctl_dir"

# Loop from start to end
for ((fam_num=start; fam_num<=end; fam_num++)); do
    fam_base="family${fam_num}"
    phy="$phy_dir/${fam_base}.cds.phy"
    ctl_file="$ctl_dir/${fam_base}.yn00.ctl"

    if [[ -f "$phy" ]]; then
        echo "▶ Generating control file for $fam_base"
        awk -v file="$phy" '{gsub("XXXXX",file); print $0}' yn00.ctl_master > "$ctl_file"
    else
        echo "Skipping $fam_base (no PHYLIP file found)"
    fi
done

echo "Control files saved in $ctl_dir for families $start to $end"
