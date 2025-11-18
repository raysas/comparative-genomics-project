## [0.0.1] - 2025-11-06 by @raysas
- connecting to github
- launching version control through ./dev/version_track.sh
- if you see this message in the changelog => test is successful, ready to go🚀

## [0.0.2] - 2025-11-15 by @raysas
- edited some folder structure to include all pipeline scripts in one place
- downloaded peptides seq for glycine max with command
- added script 0 of pipeline: extract data (peptides fasta) and summarize the protein feature in one file found in data/protein_info.csv
- added script 2 of pipeline: blast all agaisnt all and computes coverage of alignment
- added plan for scripts belonging to first chunk of pipeline: fastq -> gene families cluster
- considering filtering for longest isoform next step
- output folder hid for the moment while testing scripts
- added ncbi-blast+ as one dependency
- added pipeline plan in TODO.md
- logs folder for script outputs

## [0.0.3] - 2025-11-16 by @ayabtg
- Added 1_filter_isoforms.sh script.
- Implemented longest-isoform filtering using embedded Python.
- Generated peptides_longest.fa and protein_info_longest.csv.
- Updated the pipeline TODO list.

## [0.0.4] - 2025-11-17 by @raysas
- fixed blast
- running

## [0.0.5] - 2025-11-18 by @raysas
- fixed script 4,5,6 in pipeline_1
- refixed output directories in scripts
- got prelimnary results

