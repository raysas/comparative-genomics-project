#!/bin/bash

# - Currently hardcoded to use the "test" dataset.
# - To switch to "full" dataset, change the paths below.
input_dir="/mnt/d/documents/documents/comparative-genomics-project/comparative-genomics-project/output/ks/test/cds_pairs"
output_dir="/mnt/d/documents/documents/comparative-genomics-project/comparative-genomics-project/output/ks/test/aligned_pairs"


mkdir -p $output_dir

for file in $input_dir/*.fa; do
    base=$(basename $file .fa)
    mafft --auto $file > $output_dir/${base}_aligned.fa
done
