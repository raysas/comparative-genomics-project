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

## [0.0.6] - 2025-11-22 by @raysas
- completed yesterday's work (mainly debugging)
- changes in 5_prepare_edgelist.sh (1) and 6_cluster_families (2)
- (1) filtration and fixing output coming from script 4 (filtered blast results) to a similarity network (in form of tsv edge list)
- (1) a. made sure that the same symmetric edge (a b w = b a w) only appears once (undirected edges)
- (1) b. made sure no self loops are kept (protein self-matching, for each protein in the list)
- (1) c. removed multiple edges coming from multiple hits: for 2 proteins a and b, keep only the best hit between a and b and remove the rest (by highest edge-weight, by default, bit score)
- (2) figured an issue in MCL output: some ids are merged like KRH29797KRH39445 so read as 1 seq instead of 2
- (2) fixed it inplace and reformatted the mcl txt output into a tsv output file with the columns: geneName  family
- some fixing in input/output file names for consistency
- updated todolist for pipeline_1: officially clean pipeline
- started with small documentation, detailed one ahead
- cleaned log files from scripts

### [0.0.7] - 2025-11-22 18:00:00
- Created environment.yaml file for reproducibility.
- Created test dataset to ease debugging.
- Generated paralog pairs from clustered families.
- Produced paralog_pairs.tsv in output/ks/test/.
- Downloaded Glycine max CDS from Ensembl Plants.
- Extracted CDS sequences for paralog pairs.
- Used Bio.SeqIO to parse cds.fa and write FASTA files.
- Encountered error, changed the version of CDS, rerun.
- Created cds_pairs/ folder with individual paralog CDS files.
- Added 3_align_pairs.sh script to align CDS paralog pairs using MAFFT.
- Added inline comments to improve script readability and reproducibility.
- Added 4_calculate_ks.sh script to convert aligned FASTA to AXT and run KaKs_Calculator.
- Modified script so .axt files are saved in the same folder as aligned FASTA.
- Produced *_kaks.txt results in ks_results/ folder.
- Noted that some KaKs result files are empty due to invalid or short alignments.
- Filtered out empty *_kaks.txt files caused by invalid/short alignments.
- Merged valid Ka/Ks results, extracted Ks column and removed NA values.
- Plotted histogram of Ks distribution.
- Checked all results before implementing on full dataset.
