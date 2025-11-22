# comparative-genomics-project

[todo list for this project](TODO.md)

Pipeline 1: bash scripts located in the `pipeline_1/` folder. See [TODO.md](TODO.md) for details on pipeline steps.

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

And log files for each step can be found in the [`logs/pipeline`](logs/pipeline) directory to track the progress and debug if needed

## [0_extract_data.sh](pipeline_1/0_extract_data.sh)

> Extract protein sequences in fasta format from ENSEMBL Plants for a given species (from link) and extract relevant features into a csv file

## [1_filter_isoforms.sh](pipeline_1/1_filter_isoforms.sh)

> Filter the fasta file and features file to keep only the longest isoform per gene

## [2_blast.sh](pipeline_1/2_blast.sh)

> Perform an all-vs-all BLASTP of the protein sequences against themselves and output a tsv file with the results

## [3_compute_coverage.sh](pipeline_1/3_compute_coverage.sh)

> Compute query coverage (qcov) and subject coverage (scov) for each BLASTP hit based on sequence lengths from the features file and add these columns to the BLASTP results file

## [4_filter_pairs.sh](pipeline_1/4_filter_pairs.sh)

> Filter the BLASTP results based on user-defined thresholds for identity, qcov, and scov (default: 30, 50, 50; can also specify thresholds for bit score and e-value) and output the filtered results to a new tsv file

## [5_prepare_edgelist.sh](pipeline_1/5_prepare_edgelist.sh)

> Prepare an edgelist file from the filtered BLASTP results for clustering, allowing the user to specify which column to use as the weight (default: bit score found on the 12th column, so parameter is 12)

## [6_cluster_families.sh](pipeline_1/6_cluster_families.sh)

> Cluster the proteins into families using the MCL algorithm based on the edgelist file, allowing the user to specify MCL parameters such as inflation and whether to discard loops (self-connections)

# Requirements

* curl if not already installed
* BLAST+ installed and accessible in PATH
* MCL installed and accessible in PATH

_installation can be found in [scripts/SETUP.sh](./scripts/SETUP.sh)_