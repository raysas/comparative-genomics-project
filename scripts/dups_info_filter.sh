#!/bin/bash

# --------------------------------------------------------
# -- script to filter duplicated genes information file from the main protein info file
# -- only keeping genes that are present in the duplicated genes clusters file
# --------------------------------------------------------

# -- for low stringency clusters: (default)
# -- for high stringency clusters H:  ./scripts/dups_info_filter.sh --duplicated_genes_clusters_file output/clusters/protein_families_filtered_blast_results_id50_qcov70_scov70_evalue10e-10_wcol12_network.tsv

input_file="data/protein_info_longest.csv"
duplicated_genes_clusters_file="output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv"
# --H: output/clusters/protein_families_filtered_blast_results_id50_qcov70_scov70_evalue10e-10_wcol12_network.tsv

output_dir='output/'

usage() {
  echo "Usage: $0 [--info_file <input_file>] [--duplicated_genes_clusters_file <file>] [--output_dir <output_directory>]"
  echo "  --info_file    Input info CSV (default: $input_file)"
  echo "  --duplicated_genes_clusters_file  Duplicated genes clusters file (default: $duplicated_genes_clusters_file)"
  echo "  --output_dir   Output directory (default: $output_dir)"
}

# parse long options (only) so callers can pass --duplicated_genes_clusters_file
while [[ $# -gt 0 ]]; do
  case "$1" in
    --info_file)
      input_file="$2"; shift 2;;
    --duplicated_genes_clusters_file|--duplicated_genes_files)
      duplicated_genes_clusters_file="$2"; shift 2;;
    --output_dir)
      output_dir="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    --)
      shift; break;;
    *)
      echo "Unknown option: $1" >&2; usage; exit 1;;
  esac
done

statistics_output_dir=${output_dir}'info/'
filtration=$(basename "$duplicated_genes_clusters_file" | sed 's/protein_families_filtered_blast_results//; s/_network.tsv//')
output_file=${statistics_output_dir}'duplicated_genes_info'${filtration}'.csv'

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