#!/bin/bash

# -- script to compute distribution of duplicated genes accross chormosomes
# -- takes as input data/protein_info_longest.csv 

input_file="data/protein_info_longest.csv"
output_file='output/statitics/duplicated_genes_distribution.csv'

# -- take user input --input
usage() {
	echo "Usage: $0 [--input <input_file>]"
	echo "  --input   Path to input CSV (default: $input_file)"
}

# parse only long option --input (no short options)
while [[ $# -gt 0 ]]; do
	case "$1" in
		--input)
			input_file="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		*)
			echo "Unknown option: $1" >&2
			usage
			exit 1
			;;
	esac
done

echo "Input file: $input_file"

# (script continues...) - compute distribution using $input_file and write to $output_file