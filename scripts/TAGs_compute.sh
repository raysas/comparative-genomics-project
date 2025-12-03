#!/bin/bash

# -------------------------------------------------
# -- script to compute TAGs using the Rscript in analysis/duplicated_genes/detect_TAGs.R
# -- output files are written to output/duplication_classes/TAGs/
# -- uses default values for:
#    --dup_file: data/duplicated_genes_info_longest.csv
#    --families_file:
        # * low: data/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_evalue1e-10_wcol12_network.tsv
        # * high: data/clusters/protein_families_filtered_blast_results_id50_qcov70_scov70_evalue1e-10_wcol12_network.tsv
# -------------------------------------------------

# -- estimated run time is 4200seconds for low dataset sequentially (~70 minutes)

# > [!CAUTION]
# > nned to have some R packages installed:
# > - argparse
# > - BiocManager -> BiocManager::install("GenomicRanges")
# > - dplyr

# -- takes argument to know if low or high
mode=$1

if [[ "$mode" != "low" && "$mode" != "high" ]]; then
    echo "Usage: $0 [low|high]"
    exit 1
fi

for spacer in {0..10}
do
    echo "-- computing TAGs with spacer = $spacer for $mode stringency duplicated genes"
    (
        start=$(date +%s)
        Rscript ./scripts/TAGs_detect.R --spacer $spacer --stringency $mode
        end=$(date +%s)
        echo "-- spacer: $spacer ($mode) - elapsed time: $((end - start)) seconds"
    ) &
done

wait
echo "-- All TAGs computations for $mode stringency completed and can be found in output/duplication_classes/TAGs/$mode/"


# for spacer in {0..10}
# do
#     echo "-- computing TAGs with spacer = $spacer for low stringency duplicated genes"
#     start=$(date +%s)
#     Rscript ./scripts/TAGs_detect.R --spacer $spacer
#     end=$(date +%s)
#     echo "-- spacer: $spacer - elapsed time: $((end - start)) seconds"
# done

# for spacer in {0..10}
# do
#     echo "-- computing TAGs with spacer = $spacer for high stringency duplicated genes"
#     start=$(date +%s)
#     Rscript ./scripts/TAGs_detect.R --spacer $spacer --outdir output/duplication_classes/TAGs/ --families_file output/clusters/protein_families_filtered_blast_results_id50_qcov70_scov70_evalue1e-10_wcol12_network.tsv --dup_file  output/statistics/duplicated_genes_info_id50_qcov70_scov70_evalue1e-10_wcol12.csv
#     end=$(date +%s)
#     echo "-- spacer: $spacer - elapsed time: $((end - start)) seconds"
# done