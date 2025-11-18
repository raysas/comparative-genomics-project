# TODO

## pipeline

- [x] data extraction and protein features script [./pipeline/0_extract_data.sh](./pipeline/0_extract_data.sh)
- [x] filter isoforms script (use data/protein_info.csv to keep longest isoform per gene): [./pipeline/1_filter_isoforms.sh](./pipeline/1_filter_isoforms.sh)
- [x] all-vs-all bash script [./pipeline/2_blast.sh](./pipeline/2_blast.sh)
- [x] compute coverage script (get qcov and scov) [./pipeline/3_compute_coverage.sh](./pipeline/3_compute_coverage.sh)
- [x] filter pairs based on thresholds from blast metrics [pipeline_1/4_filter_pairs.sh](./pipeline/4_filter_pairs.sh) : add option to filter based on different metrics (identity, qcov, scov, bit score, evalue) - user defined thresholds
- [x] get edgelist file (putative homologs with bit score) script | filter options based on different metrics (user defines) [./pipeline/5_prepare_edgelist.sh](./pipeline/5_prepare_edgelist.sh)
- [x] clustering script (check which algos - ones that match the galaxy ones) | get gene families [./pipeline/6_cluster_families.sh](./pipeline/6_cluster_families.sh)
- [ ] get another format for clusters: 2 columns for gene and family ID (7th script)

from the end file later on perform filtering, analysis, plotting, stats etc

- [ ] make a nice diagram of the pipeline steps
- [ ] finalize bash version of the pipeline
- [ ] create a docker container for pipeline dependencies and link it to dockerhub
- [ ] modularize it for nextflow (+ base it on nextflow image)