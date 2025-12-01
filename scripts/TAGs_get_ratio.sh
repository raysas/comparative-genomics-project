#!/bin/bash

# ------------------------------------------------------
# -- this scripts collectes inforamtion about TAGs retrieved with different spacers on:
#   - a=number of TAG genes retrieved
#   - b=number of TAGs arrays retrieved
#   - c=ratio of TAG genes / total number of genes
# -- in an attempt to replicate table in Shoja & Zhang, 2006
# ------------------------------------------------------

output_dir='output/statistics/'
output_file=${output_dir}'TAGs_ratios.tsv'

cluster_dir='output/duplication_classes/TAGs/'
input_file='data/protein_info_longest.csv' # --> to get total number of genes

total_genes=$(tail -n +2 $input_file | wc -l)
echo -e "spacer\tn_TAG_genes\tn_TAG_arrays\tratio_TAG_genes" > $output_file

for file in ${cluster_dir}TAGs_*.tsv
do
    echo "-- processing file: $file"
    spacer=$(basename $file | sed 's/TAGs_\(.*\).tsv/\1/')
    # -- number of tag genes is how much rows in 6th column>=1
    n_TAG_genes=$(awk -F'\t' 'NR>1 {if ($6 >=1) count++} END {print count+0}' $file)
    # -- number of tag arrays is how much unique 5th and 6th column combinations
    n_TAG_arrays=$(awk -F'\t' 'NR>1 {if ($6 >=1) arr[$5]++} END {print length(arr)+0}' $file)
    ratio_TAG_genes=$(echo "scale=5; $n_TAG_genes / $total_genes" | bc)
    echo -e "${spacer}\t${n_TAG_genes}\t${n_TAG_arrays}\t${ratio_TAG_genes}" >> $output_file
done

echo "-- results written to: $output_file"