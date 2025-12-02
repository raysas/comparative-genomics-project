#!/bin/bash

# -------------------------------------------------
# -- script to compute TAGs using the Rscript in analysis/duplicated_genes/detect_TAGs.R
# -- output files are written to output/duplication_classes/TAGs/
# -- uses default values for:
#    --dup_file: data/duplicated_genes_info_longest.csv
#    --families_file: data/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv
# -------------------------------------------------

# -- estimated run time is 4200seconds

# > [!CAUTION]
# > nned to have some R packages installed:
# > - argparse
# > - BiocManager -> BiocManager::install("GenomicRanges")
# > - dplyr


for spacer in {0..10}
do
    echo "-- computing TAGs with spacer = $spacer"
    start=$(date +%s)
    Rscript ./scripts/detect_TAGs.R --spacer $spacer --outfile output/duplication_classes/TAGs/TAGs_${spacer}.tsv
    end=$(date +%s)
    echo "-- spacer: $spacer - elapsed time: $((end - start)) seconds"
done