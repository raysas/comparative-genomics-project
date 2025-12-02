#!/bin/bash

# --------------------------------------------------------
# -- script to filter duplicated genes information file from the main protein info file
# -- only keeping genes that are present in the duplicated genes clusters file
# --------------------------------------------------------


input_file="data/protein_info_longest.csv"
duplicated_genes_clusters_file="output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv"

output_dir='output/'

while getopts ":info_file:duplicated_genes_files:output_dir:" opt; do
  case $opt in
    info_file)
      input_file="$OPTARG"
      ;;
    duplicated_genes_files)
      duplicated_genes_clusters_file="$OPTARG"
      ;;
    output_dir)
      output_dir="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      ;;
  esac
done

statistics_output_dir=${output_dir}'statistics/'
output_file=${statistics_output_dir}'duplicated_genes_info.csv'

if [ ! -d "$statistics_output_dir" ]; then
    mkdir -p "$statistics_output_dir"
fi

echo " -- input file: "$input_file
echo " -- duplicated genes clusters file: "$duplicated_genes_clusters_file

# -- first copying the input info file in statitics output dir for reference
cp "$input_file" "${statistics_output_dir}$(basename "$input_file")"

# -- filtering for genes present in the clusters file
header=$(head -n 1 $input_file)
echo "$header" > $output_file
awk -F'\t' 'NR==FNR { ids[$1]; next }
            {
                split($0, a, ",");
                if (a[1] in ids) print $0
            }' \
    $duplicated_genes_clusters_file \
    $input_file >> $output_file
echo " -- filtered duplicated genes info written to: "$output_file