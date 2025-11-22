# comparative-genomics-project: pipeline_1

[todo list for this project](TODO.md)

Pipeline 1: bash scripts located in the `pipeline_1/` folder. See [TODO.md](TODO.md) for details on pipeline steps (in branch `pipeline_1`)

# Documentation

This pipeline is made of 6 main steps +1 for extracting peptides if not already provided:

```
pipeline_1/
├── 0_extract_data.sh
├── 1_filter_isoforms.sh
├── 2_blast.sh
├── 3_compute_coverage.sh
├── 4_filter_pairs.sh
├── 5_prepare_edgelist.sh
└── 6_cluster_families.sh
```

Will mention for each file:  
* Input files/arguments
* Output files/directory
* Main parameters
* Short description of the step and problems that were solved

To extract _Glycine max_ protein sequences and features from ENSEMBL Plants, can run:
```bash
./pipeline_1/0_extract_data.sh      # -- default param for Glycine max
./pipeline_1/1_filter_isoforms.sh   # -- to filer longest isoforms only
```
Will return the following directory:
```
data/
├── peptides.fa.gz
├── peptides.fa
├── peptides_longest.fa         # -- longest isoforms only 
├── protein_info.csv            # -- MAIN FEATURES FILE: ids, chromosome, start, end, strand, description seq length (concat of metadata+seqlength)
├── protein_info_longest.csv    # -- kept all protein features for longest isoforms only
├── protein_metadata.csv        # -- metadata extracted from fasta headers (chr number, start, end, strand, and ids)
└── seq_lengths.csv             # -- sequence lengths for all proteins
```

> [!TIP]
> (later fix)  
> better to specify the default directory where this will be downlaoded as `data/Glycine_max/` for better organization when adding more species, not changing this now to avoid breaking the pipeline (also then make `output/Glycine_max/` for outputs instead of just `output/`)

To run the full pipeline with default values can do the following:
```bash
./pipeline_1/2_blast.sh             # -- all-vs-all blastp, if output/blst_output exists, will skip
./pipeline_1/3_compute_coverage.sh  # -- compute qcov and scov from seq length prev computed in script 0 and add to blast output
./pipeline_1/4_filter_pairs.sh      # -- filter based on identity, qcov, scov thresholds (default: 30, 50, 50)
./pipeline_1/5_prepare_edgelist.sh  # -- prepare edgelist file from filtered blast results, can specify weight column
./pipeline_1/6_cluster_families.sh  # -- cluster using mcl, can specify mcl paramemeters
```

> [!NOTE]
> Each script has arguments to be modified, will explain them down below (e.g., thresholds for filtering, input/output files etc)

Will get this as output:
```
output/
├── blast_filtered          # -- script 4 => filtered blast results based on thresholds
│   └── filtered_blast_results_id30_qcov50_scov50.tsv
├── blast_output/           # -- script 2 => all-vs-all blast results
│   ├── blast_results.tsv
│   ├── blast_results_with_coverage.tsv   # -- script 3 => blast results with coverage columns added
│   ├── peptide_db.pdb
│   ├── peptide_db.phr
│   ├── peptide_db.pin
│   ├── peptide_db.pot
│   ├── peptide_db.psq
│   ├── peptide_db.ptf
│   └── peptide_db.pto
├── clusters                # -- script 6 => clustering results (protein families)
│   ├── protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12.tsv # -- reformatted
│   └── protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12.txt # -- mcl output format
└── similarity_edgelists    # -- script 5 => edgelist (similarity network) files for clustering
    └── filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv
```

And log files for each step can be found in the [`logs/pipeline`](logs/pipeline) directory to track the progress and debug if needed.
For the structure of the script there will be mainly 3 well-commetned chunks to facilitate understanding and modification:
1. documentation and usage (comments)
2. preparing variables, getting arguments, setting up logging
3. main script steps

## [0_extract_data.sh](pipeline_1/0_extract_data.sh)

> Extract protein sequences in fasta format from ENSEMBL Plants for a given species (from link) and extract relevant features into a csv file

Arguments: 
```
Usage: ./pipeline_1/0_extract_data.sh [-l download_link] [-o output_directory] [-i FASTA_FILE]
  -l download link for peptide fasta file (default: http://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.0.pep.all.fa.gz)
  -o output directory to save  data files in (default: data)
  -i fasta file name to be saved (default: peptides.fa)
```
<!-- # --------------------------------------------------------------------
# -- What this script does:
#    1) downloads peptide data from Ensembl Plants
#    2) extracts protein information from fasta headers
#       and saves it in a tabular format: 2 outputs
#         a- the fasta file with peptides
#         b- a tsv file with extracted protein info
# -- temporary files created during processing will be stored in the output directory
#    1. metadata.csv
#    2. seq_lengths.csv
# -- Usage:
#    bash ./pipeline_1/0_extract_data.sh [-l download_link] [-o output_directory] [-i FASTA_FILE]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/0_extract_data.sh -l "http://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.0.pep.all.fa.gz" -o "data" -i "peptides.fa"
# --------------------------------------------------------------------
 -->

**Input**:
* LINK: download link for peptide fasta file (default: ftp://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.0.pep.all.fa.gz)
* OUTPUT_DIR: output directory to save data files in (default: data/)
* FASTA_FILE: fasta file name to be saved (default: peptides.fa)

**Output** - will produce several files in the output directory:
* peptides.fa.gz : downloaded compressed fasta file  
* peptides.fa : decompressed fasta file
* protein_metadata.csv : extracted protein metadata from fasta headers
* seq_lengths.csv : sequence lengths for each protein
* protein_info.csv : main features file with all info (metadata + sequence lengths)

*protein_info.csv format:*

| peptide_id | gene_id | transcript_id | genome | chromosome | start_pos | end_pos | strand | description | length |
|------------|---------|---------------|--------|------------|-----------|---------|--------|-------------|--------|


## [1_filter_isoforms.sh](pipeline_1/1_filter_isoforms.sh)

> Filter the fasta file and features file to keep only the longest isoform per gene
>
<!-- > # --------------------------------------------------------------------
# -- What this script does:
#    keeps only the longest isoform per gene from the initial peptide fasta file
#    and creates:
#      1) a filtered fasta file (e.g., peptides_longest.fa)
#      2) a filtered protein info file (e.g., protein_info_longest.csv)
# -- Usage:
#    bash ./pipeline_1/1_filter_isoforms.sh [-f FASTA_FILE] [-i FEATURE_FILE]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/1_filter_isoforms.sh -f "data/peptides.fa" -i "data/protein_info.csv"
# --------------------------------------------------------------------
###########################################################################

# ---------------------------------------------------------------------
# prepare variables, get arguments, set up logging
# ---------------------------------------------------------------------
 -->

Arguments:
```
$ ./pipeline_1/1_filter_isoforms.sh -h  
Usage: ./pipeline_1/1_filter_isoforms.sh [-f FASTA_FILE] [-i FEATURE_FILE]
  -f    Input peptide FASTA file (default: data/peptides.fa)
  -i    Input protein info CSV file (default: data/protein_info.csv)
  -h    Show this help message
```
**Input**:
* FASTA_FILE: Input peptide FASTA file (default: data/peptides.fa)  
* FEATURE_FILE: Input protein info CSV file (default: data/protein_info.csv)

**Output** - will produce the following files in the same directory as the input files:
* peptides_longest.fa : filtered fasta file with only longest isoforms per gene
* protein_info_longest.csv : filtered protein info CSV file with only longest isoforms per gene

## [2_blast.sh](pipeline_1/2_blast.sh)

> Perform an all-vs-all BLASTP of the protein sequences against themselves and output a tsv file with the results

Arguments:
```
$./pipeline_1/2_blast.sh -h
Usage: ./pipeline_1/2_blast.sh [-i input_file] [-o output_directory] [-f] [-h]
  -i    Input file
  -o    Output directory
  -f    Force creation of a new blast_output directory (do not reuse existing one - wont delete previous one will only rename it with prefix old_ and a number suffix)
  -h    Show this help message
```


**Input**:
* INPUT_FILE: Input peptide FASTA file (default: data/peptides_longest.fa)  
* OUTPUT_DIR: Output directory for BLAST results (default: output/blast_output)
* FORCE_NEW: Force creation of a new blast_output directory (do not reuse existing one - wont delete previous one will only rename it with prefix old_ and a number suffix) (default: false)

> [!NOTE]
> Did not make parameters for BLAST itself (e.g., evalue, max target seqs etc) to keep it simple, can modify the script directly if needed

**Output** - will produce the following files in the output directory:
* blast_results.tsv : BLASTP results in tsv format (tab-separated values)
* BLAST database files (various extensions: .pdb, .phr, .pin, .pot, .psq, .ptf, .pto)

*blast_results.tsv output, notice no header*

| KRH43692 | KRH43692 | 100.000 | 396 | 0  | 0 | 1 | 396 | 1 | 396 | 0.0 | 827 |
|--|---------|----------|---------|--------|----|---|---|-----|---|-----|------|-----|
| KRH43692 | KRH13762 | 94.444  | 396 | 22 | 0 | 1 | 396 | 1 | 396 | 0.0 | 780 |
| KRH43692 | KRH13763 | 89.421  | 397 | 41 | 1 | 1 | 396 | 1 | 397 | 0.0 | 737 |
| KRH43692 | KRH37933 | 77.805  | 410 | 77 | 3 | 1 | 396 | 1 | 410 | 0.0 | 650 |
| KRH43692 | KRH38261 | 69.067  | 375 | 112| 1 | 22| 396 | 29| 399 | 0.0 | 538 |


## [3_compute_coverage.sh](pipeline_1/3_compute_coverage.sh)

> Compute query coverage (qcov) and subject coverage (scov) for each BLASTP hit based on sequence lengths from the features file and add these columns to the BLASTP results file


Arguments:
```
$ ./pipeline_1/3_compute_coverage.sh -h
Usage: ./pipeline_1/3_compute_coverage.sh [-i input_file] [-o coverage_output_filename] [-p protein_info_file]
  -i    Input BLAST results file
  -o    Output file name for coverage results
  -p    Protein info file
  -h    Show this help message
```

**Input**:
* INPUT_FILE: Input BLAST results file (default: output/blast_output/blast_results.tsv)  
* COVERAGE_OUTPUT_FILE: Output file name for coverage results (default: output/blast_output/blast_results_with_coverage.tsv)  
* PROTEIN_INFO_FILE: Protein info file (default: data/protein_info_longest.csv)

**Output** - will produce the following file:
* blast_results_with_coverage.tsv : BLASTP results file with added columns for query coverage (qcov) and subject coverage (scov)

*blast_results_with_coverage.tsv output, with header:*

| qseqid  | sseqid   | pident  | length | mismatch | gapopen | qstart | qend | sstart | send | evalue     | bitscore | qlength | slength | qcov    | scov    |
|---------|----------|---------|--------|----------|---------|--------|------|--------|------|------------|----------|---------|---------|---------|---------|
| KRG58484 | KRG58484 | 100.000 | 158    | 0        | 0       | 1      | 158  | 1      | 158  | 1.12e-118  | 332      | 158     | 158     | 100     | 100     |
| KRG58484 | KRG95435 | 49.275  | 69     | 25       | 5       | 179    | 242  | 34     | 97   | 4.00e-06   | 46.2     | 158     | 255     | 40.5063 | 25.098  |
| KRG58484 | KRH38054 | 81.731  | 104    | 13       | 4       | 194    | 297  | 1      | 98   | 1.07e-30   | 113      | 158     | 307     | 65.8228 | 31.9218 |
| KRG58484 | KRH67106 | 48.387  | 62     | 23       | 4       | 184    | 240  | 35     | 92   | 2.04e-06   | 47.0     | 158     | 259     | 36.0759 | 22.3938 |


## [4_filter_pairs.sh](pipeline_1/4_filter_pairs.sh)

> Filter the BLASTP results based on user-defined thresholds for identity, qcov, and scov (default: 30, 50, 50; can also specify thresholds for bit score and e-value) and output the filtered results to a new tsv file

Arguments:
```
$ ./pipeline_1/4_filter_pairs.sh -h
Usage: ./pipeline_1/4_filter_pairs.sh [-i input_file] [-o output_dir] [-id identity_threshold] [-qcov query_coverage_threshold] [-scov subject_coverage_threshold] [-bit bit_score_threshold] [-evalue e_value_threshold]
  -i        Input BLAST results file with coverage
  -o        Output directory for filtered results
  -id       Identity threshold (default: 30)
  -qcov     Query coverage threshold (default: 50)
  -scov     Subject coverage threshold (default: 50)
  -bit      Bit score threshold (default: no filtering)
  -evalue   E-value threshold (default: no filtering)
```

**Input**:
* INPUT_FILE: Input BLAST results file with coverage (default: output/blast_output/blast_results_with_coverage.tsv)  
* OUTPUT_DIR: Output directory for filtered results (default: output/blast_filtered)
* ID_THRESHOLD: Identity threshold (default: 30)  
* Q_COV_THRESHOLD: Query coverage threshold (default: 50)   
* S_COV_THRESHOLD: Subject coverage threshold (default: 50)
* BIT_SCORE_THRESHOLD: Bit score threshold (default: no filtering)  
* E_VALUE_THRESHOLD: E-value threshold (default: no filtering)

**Output** - will produce the following file in the output directory:
* filtered_blast_results_id{ID_THRESHOLD}_qcov{Q_COV_THRESHOLD}_scov{S_COV_THRESHOLD}.tsv : filtered BLASTP results based on the specified thresholds

## [5_prepare_edgelist.sh](pipeline_1/5_prepare_edgelist.sh)

> Prepare an edgelist file from the filtered BLASTP results for clustering, allowing the user to specify which column to use as the weight (default: bit score found on the 12th column, so parameter is 12)

Arguments:
```
$ ./pipeline_1/5_prepare_edgelist.sh -h
Usage: ./pipeline_1/5_prepare_edgelist.sh [-i input_file] [-o OUTPUT_DIR] [-w weight_column_index]
  -i    Input filtered BLAST results file
  -o    Output edgelist file
  -w    Column index for weight in BLAST output
  -h    Show this help message
```


**Input**:
* INPUT_FILE: Input filtered BLAST results file (default: output/blast_filtered/filtered_blast_results_id30_qcov50_scov50.tsv)  
* OUTPUT_DIR: Output edgelist file (default: output/similarity_edgelists)  
* WEIGHT_COLUMN_INDEX: Column index for weight in BLAST output (default: 12, which corresponds to bit score)

**Output** - will produce the following file in the output directory:
* filtered_blast_results_id{ID_THRESHOLD}_qcov{Q_COV_THRESHOLD}_scov{S_COV_THRESHOLD}_wcol{WEIGHT_COLUMN_INDEX}_network.tsv : edgelist file for clustering

*filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv  format: protein1 protein2 weight*

| KRG58484 | KRG95435 | 46.2  |
|--|---------|-------|
| KRG58484 | KRH38054 | 113   |
| KRG58484 | KRH67106 | 47.0  |
| KRG61172 | KRG97777 | 29.3  |
| KRG62078 | KRH33086 | 154   |
| KRG62828 | KRG93564 | 82.4  |
| KRG62828 | KRG93574 | 266   |



## [6_cluster_families.sh](pipeline_1/6_cluster_families.sh)

> Cluster the proteins into families using the MCL algorithm based on the edgelist file, allowing the user to specify MCL parameters such as inflation and whether to discard loops (self-connections)

Arguments:
```
$ ./pipeline_1/6_cluster_families.sh -h
Usage: ./pipeline_1/6_cluster_families.sh [-i input_file] [-o output_dir]
  -i    Input edgelist file for clustering
  -o    Output file for protein families
  -h    Show this help message
  -c <num>           increase loop-weights <num>-fold (default: 1.0)
  -discard_loops y|n   discard self-loops (default: y)
  -I <num>          MCL inflation parameter (default: 2.0)
  -P <num>          MCL (inverted) rigid pruning threshold (cf -z) (default: 4000)
  -S <num>          MCL select down to <int> entries if needed
  -R <num>          MCL recover to maximally <int> entries if needed
  -pct <num>        MCL try recovery if mass is less than <pct>
for more info see mcl --help
```


> [!CAUTION]
> there was an issue with merged ids (e.g., KRH29797KRH39445) in the output of MCL, fixed now through a python script using regex library inside the bash script, i.e. both outputs will have the fixed ids

**Input**
* INPUT_FILE: Input edgelist file for clustering (default: output/similarity_edgelists/filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv)  
* OUTPUT_DIR: Output file for protein families (default: output/clusters)
* MCL parameters: can specify various MCL parameters such as inflation (-I), discard loops (-discard_loops), loop weight (-c), pruning threshold (-P), select (-S), recover (-R), and pct (-pct). Default values are provided in the script.

**Output** - will produce the following files in the output directory:
* protein_families_{input_filename}.txt : MCL output format with each line representing a family with space-separated gene IDs (fixed merged ids issue)
* protein_families_{input_filename}.tsv : TSV format with two columns: geneName and familyID

*`.txt` (protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.txt) file format*

| KRH43692 KRH13762 KRH13763 KRH37933 KRH38261 |
|-----------------------------------------------|
| KRG58484 KRG95435 KRH38054 KRH67106           |
| KRG61172 KRG97777                             |

*`.tsv` (protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv) file format*

|geneName| family|
|--------|-------|
|KRH43692| 1     |
|KRH13762| 1     |
|KRH13763| 1     |
|KRH37933| 1     |
|KRH38261| 1     |

> [!NOTE]
> Also worth noting, the script now removes clusters with only 1 member (not interested in singleton clusters - no gene duplicates information), so the final output files will only contain clusters with 2 or more members (minimum size of a gene family)

# Requirements

* curl if not already installed
* BLAST+ installed and accessible in PATH
* MCL installed and accessible in PATH

_installation can be found in [scripts/SETUP.sh](./scripts/SETUP.sh)_