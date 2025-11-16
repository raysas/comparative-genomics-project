# TODO

## pipeline

- [x] data extraction and protein features script [./pipeline/0_extract_data.sh](./pipeline/0_extract_data.sh)
- [x] filter isoforms script (use data/protein_info.csv to keep longest isoform per gene): [./pipeline/1_filter_isoforms.sh](./pipeline/1_filter_isoforms.sh)
- [x] all-vs-all bash script [./pipeline/2_blast.sh](./pipeline/2_blast.sh)
- [ ] get edgelist file (putative homologs with bit score) script | filter options based on different metrics (user defines) [./pipeline/3_filter_pairs.sh](./pipeline/3_filter_pairs.sh)
- [ ] clustering script (check which algos - ones that match the galaxy ones) | get gene families [./pipeline/4_cluster_families.sh](./pipeline/4_cluster_families.sh)
- [ ] get gene pairs within cluster/fam script (same thing as above ig)
- [ ] get seq of pair of genes scritp (prep for next step, maybe merge with)
- [ ] compute Ks wirth PAML script (get file with Ks values)

from the end file later on perform filtering, analysis, plotting, stats etc

- [ ] make a nice diagram of the pipeline steps
- [ ] finalize bash version of the pipeline
- [ ] create a docker container for pipeline dependencies and link it to dockerhub
- [ ] modularize it for nextflow (+ base it on nextflow image)