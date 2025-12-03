#!/bin/bash

# -- this script is sourced by scripts/dups_singletons_extraction.r
# -- it takes as input:
#    - duplication families cluster file path
#    - output dir path
# OR
#    - stringency level: low or high, then will take care of setting the file paths accordingly
# -- it extracts singletons from dups directly
# -- also saves singletons info file in output/info/

# -- get args using getopts

dup_file=""
outdir=""
stringency=""

while getopts d:o:s:h flag
do
    case "${flag}" in
        d) dup_file=${OPTARG};;
        o) outdir=${OPTARG};;
        s) stringency=${OPTARG};;
        h) echo "Usage: $0 [-d duplication_families_file] [-o output_directory] [-s stringency_level (low|high)]"
           exit 0;;
        *) echo "Invalid option: -${flag}" >&2
           exit 1;;
    esac
done

# need to have either dup_file and outdir or stringency
if [[ -z "$dup_file" || -z "$outdir" ]]; then
    if [[ -z "$stringency" ]]; then
        echo "Error: Either provide both --dup_file and --outdir, or provide --stringency (low or high)."
        exit 1
    fi
fi

if [[ -n "$stringency" ]]; then
    echo "-- stringency option provided: $stringency, overriding any provided file paths"
    if [[ "$stringency" == "low" ]]; then
        echo "-- using low stringency duplicated genes and families, overriding any provided file paths: id30, cov50, evalue1e-10"
        dup_file="output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_evalue1e-10_wcol12_network.tsv"
        outdir="output/gene_lists/singletons/"
    elif [[ "$stringency" == "high" ]]; then
        echo "-- using high stringency duplicated genes and families, overriding any provided file paths: filter id50, cov70, evalue1e-10"
        dup_file="output/clusters/protein_families_filtered_blast_results_id50_qcov70_scov70_evalue1e-10_wcol12_network.tsv"
        outdir="output/gene_lists/singletons/"
    else
        echo "Error: Invalid stringency level. Use 'low' or 'high'."
        exit 1
    fi
fi

echo "-- duplication families file: $dup_file"
echo "-- output directory: $outdir"
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir"
fi

# -- extract singletons from total
prot_info="data/protein_info_longest.csv"
# peptide_id,gene_id,transcript_id,genome,chromosome,start_pos,end

# -- 1. gene_lists/ output

peptide_ids_in_dups=$(awk -F'\t' '{for(i=1;i<=NF;i++) print $i}' "$dup_file" | sort | uniq)
echo "-- number of unique peptide IDs in duplicated families: $(echo "$peptide_ids_in_dups" | wc -l)"
# total_peptide_ids=$(tail -n +2 "$prot_info" | cut -d',' -f1 | sort | uniq)
# singletons=$(comm -13 <(echo "$peptide_ids_in_dups") <(echo "$total_peptide_ids"))
# echo "$singletons" > "${outdir}singletons_${stringency}.txt"
# cat $prot_info | grep -v  awk -F',' '{print $1}' > "${outdir}singletons_${stringency}.txt"

echo "-- singleton peptide IDs extracted and saved to ${outdir}singletons_${stringency}.txt"

# 2.
output_file="output/info/singletons_genes_info_${stringency}.csv"
header=$(head -n 1 $prot_info)
echo "$header" > $output_file

awk -F'\t' '
    FNR==NR {
        for(i=1;i<=NF;i++) dups[$i]=1
        next
    }
    FNR>1 { 
        if(!dups[$1]) print $1
    }
' "$dup_file" "$prot_info" >> $output_file
echo "-- singleton gene info saved to $output_file"