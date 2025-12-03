#!/bin/bash


# ./scripts/dups_get_list.sh # -- for low
# ./scripts/dups_get_list.sh high # -- for high

stringency=low
if [ "$1" == "high" ]; then
    stringency="high"
fi

input_file="output/info/duplicated_genes_info_${stringency}.csv"
output_file="output/gene_lists/duplicated_genes/duplicated_genes_${stringency}.txt"

echo "-- extracting duplicated genes list from: $input_file"
if [ ! -d "output/gene_lists/duplicated_genes/" ]; then
    mkdir -p "output/gene_lists/duplicated_genes/"
fi

tail -n +2 $input_file | awk -F ',' '{print $1}' > $output_file
echo "-- duplicated genes list saved to: $output_file"
