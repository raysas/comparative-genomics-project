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

