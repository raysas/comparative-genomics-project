# Gene Duplication and Transposable Elements in *Glycine max*

![Domain](https://img.shields.io/badge/comparative%20genomics-plant-green)
![Pipeline](https://img.shields.io/badge/pipeline-passing-brightgreen)
![Status](https://img.shields.io/badge/status-course_project-orange)  
![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white)
![Perl](https://img.shields.io/badge/Perl-39457E?logo=perl&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-4EAA25?logo=gnubash&logoColor=white)


## Overview
Computational genomics analysis of gene duplication mechanisms and transposable element (TE) dynamics in the *Glycine max* genome.  
The project integrates large-scale sequence clustering, Ks-based duplication dating, statistical testing, and comparative genomics to study the impact of whole-genome duplication, tandem duplication, and TE activity on soybean genome evolution.

---

## Data
- *Glycine max* genome and annotations (EnsemblPlants v2.0)
- Protein and CDS FASTA files
- GO Slim functional annotations (PANTHER v19.0)
- Transposable element annotations (APTEdb)

---

## Pipeline
The analysis is implemented as a **fully scripted, modular pipeline** to generate gene families and Ks values.  
Will provide reproducible instructions for running the pipeline on new datasets.



## Analyses
- Protein family inference using BLASTP and MCL
- Ks estimation via protein-guided codon alignments (MAFFT, PAL2NAL, PAML)
- Tandemly arrayed gene (TAG) detection and orientation analysis
- Functional enrichment analysis (Fisher’s exact test, FDR correction)
- Transposable element abundance, coverage, and nesting analyses
- Gene–TE spatial proximity analysis
- Statistical testing and visualization in R and Python

---

## Datasets
Two duplication datasets were generated using different stringency thresholds:

- **Low stringency** (`low` / `id30_cov50_evalue1e-10`)  
  - 30% sequence identity  
  - 50% query and subject coverage  
  - e-value = 1e-10  

- **High stringency** (`high` / `id50_cov70_evalue1e-10`)  
  - 50% sequence identity  
  - 70% query and subject coverage  
  - e-value = 1e-10  

Both datasets were analyzed independently to assess robustness and parameter sensitivity.