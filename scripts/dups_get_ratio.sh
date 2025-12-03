#!/bin/bash

# ------------------------------------------------------
# -- script aims to get # of duplicated genes / total genes ratio
# -- takes all files in output/clusters/ to check the ratio using different thresholds for filtering
# ------------------------------------------------------


output_dir='output/statistics/'
output_file=${output_dir}'duplication_ratios.tsv'

cluster_dir='output/clusters/'
input_file='data/protein_info_longest.csv'

total_genes=$(($(wc -l < $input_file) - 1))  #1 for header

echo -e "cluster_file\tduplicated_genes\ttotal_genes\tduplication_ratio\tlargest_family_size" > $output_file

for cluster_file in ${cluster_dir}*.tsv
do
    echo "-- processing file: $cluster_file for duplicated genes ratio"
    duplicated_genes=$(($(awk '{print $1}' $cluster_file | sort | uniq | wc -l)-1))
    duplication_ratio=$(echo "scale=4; $duplicated_genes / $total_genes" | bc)

    echo "-- processing file: ${cluster_file%.tsv}.txt for largest family size"
    text_cluster_file=${cluster_file%.tsv}.txt
    largest_family_size=$(head -n 1 $text_cluster_file | awk '{print NF}')
    echo ""

    echo -e "$(basename $cluster_file)\t$duplicated_genes\t$total_genes\t$duplication_ratio\t$largest_family_size" >> $output_file
done

echo "-- duplication ratios saved to $output_file:"

cat $output_file