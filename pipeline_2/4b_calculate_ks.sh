#!/bin/bash

# - Currently hardcoded to use the "test" dataset.
# - To switch to "full" dataset, change the paths below.
input_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks_b/test/aligned_pairs"
output_dir="/mnt/d/documents/Nhi/comparative-genomics-project/output/ks_b/test/ks_results"
mkdir -p "$output_dir"

# -- Loop through all aligned FASTA files
for file in "$input_dir"/*_aligned.fa; do
    base=$(basename "$file" _aligned.fa)

    # -- Save .axt in the same folder as the FASTA file
    axt_file="${input_dir}/${base}.axt"

    # -- Save Ka/Ks results in ks_results folder
    result_file="${output_dir}/${base}_kaks.txt"

    # -- Convert aligned FASTA to AXT format (required by KaKs_Calculator)
    awk '/^>/{if (NR==1) {printf "%s\t", substr($0,2)} else {printf "\n%s\t", substr($0,2); next} next} {printf "%s", $0} END{print ""}' "$file" > "$axt_file"

    # -- Run KaKs_Calculator with NG model
    KaKs_Calculator -i "$axt_file" -o "$result_file" -m NG
done
