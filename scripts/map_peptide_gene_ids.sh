#!/bin/bash

# ------------------------------------------------------------
# script to map identifiers between peptide and gene ids
# take a protein info file as reference (data/protein_info_longest.csv)
# takes (2 options):
# 1) --input_file : path to the gene list file (txt file like the ones in gene_lists)
# 2) --input_dir : path to a directory containing multiple gene list files (txt files like the ones in gene_lists)
# 3) --peptide : just a peptide id(s) as argument (sysout output)
# 4) --outdir : output directory (default: output/gene_lists/mapped_peptide_gene_ids/)
# [!IMPORTANT] output is only gene ids no peptide ids 
# ------------------------------------------------------------

prot_info="data/protein_info_longest.csv"
# peptide_id,gene_id,transcript_id,genome,chromosome,start_pos,end_pos,strand,description,length

input_file=""
input_dir=""
peptide_ids=()
output_dir="output/gene_lists/mapped_peptide_gene_ids/"

usage() {
  echo "Usage: $0 [--input_file <file>] [--input_dir <directory>] [--peptide <peptide_id1 peptide_id2 ...>] [--outdir <output_directory>]"
  echo "  --input_file   Path to a gene list file (txt)"
  echo "  --input_dir    Path to a directory containing multiple gene list files (txt)"
  echo "  --peptide      One or more peptide IDs to map (space-separated)"
  echo "  --outdir       Output directory (default: $output_dir)"
}

# parse long options
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input_file)
      input_file="$2"; shift 2;;
    --input_dir)
      input_dir="$2"; shift 2;;
    --peptide)
      shift
      while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
        peptide_ids+=("$1")
        shift
      done
      ;;
    --outdir)
      output_dir="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    --)
      shift; break;;
    *)
      echo "Unknown option: $1" >&2; usage; exit 1;;
  esac
done

if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# -- option 1
if [ -n "$input_file" ]; then
    echo "-- mapping peptide to gene IDs for file: $input_file"
    base_name=$(basename "$input_file" .txt)
    output_file="${output_dir}${base_name}_mapped.txt"

    # prot_info is CSV with columns: peptide_id,gene_id,...
    # Use a flexible field separator to handle comma or tab
    awk -F'[\t,]' 'NR==FNR { map[$1]=$2; next } { key=$1; if (key in map) print map[key]; else print key"\tNA" }' "$prot_info" "$input_file" > "$output_file"
    echo "Output written to: $output_file"
fi

# -- option 2
if [ -n "$input_dir" ]; then
    echo "-- mapping peptide to gene IDs for all files in directory: $input_dir"
    output_dir="${output_dir}$(input_dir##*/)/"
    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
    fi
    for file in "$input_dir"/*.txt; do
        echo "-- processing file: $file"
        base_name=$(basename "$file" .txt)
        output_file="${output_dir}${base_name}_mapped.txt"

        awk -F'[\t,]' 'NR==FNR { map[$1]=$2; next } { key=$1; if (key in map) print map[key]; else print key"\tNA" }' "$prot_info" "$file" > "$output_file"
        echo "Output written to: $output_file"
    done
fi

# -- option 3
if [ ${#peptide_ids[@]} -gt 0 ]; then
    echo "-- mapping provided peptide IDs to gene IDs"
    # build map from prot_info then lookup each provided peptide id
    declare -A _map
    while IFS=',' read -r pep gid _rest; do
      _map["$pep"]="$gid"
    done < "$prot_info"
    for pid in "${peptide_ids[@]}"; do
      gene_id="${_map[$pid]}"
      if [ -n "$gene_id" ]; then
        echo -e "$gene_id"
      else
        echo -e "$pid\tNA"
      fi
    done
fi