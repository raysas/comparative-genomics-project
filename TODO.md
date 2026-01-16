# TODO - Comparative Genomics Project

## Major Themes

- [ ] pipeline 1: from raw data to duplicated genes
- [ ] pipeline 2: from duplicated genes to Ks computation
- [ ] duplicated genes classification (TAGs, WGD...)
- [ ] comparative analysis of singletons and duplicated genes families
- [ ] analysis of Ks distributions and age of duplications
- [ ] TEs annotation and analysis
- [ ] Computation of TE coverage


### Extensions

- [ ] Transcriptomics data integration into the project: check [this dataset](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-4279?query=glycine%20max) from arrayExpress for RNA-seq data of Glycine max - accession:E-MTAB-4279
- [ ] Comparative orthologs analysis: based on fully reproducible one-liner pipeline integration
- [ ] Functional analysis of duplicated genes: at least GO enrichment (could be on panther or using R packages like clusterProfiler, it should be following the same ORA method based on statistical tests)
- [ ] Localization and orientation of TAGs (also can do some functional/ppi types of analysis to compare)
- [ ] retrieval of PPIs from STRING database from a list of genes

## Project Overview
**Organism:** *Glycine max* (Soybean)  
**Data Source:** Ensembl Plants Release 41 - http://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/  
**Goal:** Analyze duplicated genes, compute evolutionary metrics (Ks values), and investigate transposable elements and tandemly arrayed genes

---

## Core Pipeline Tasks

### 1. Pipeline Part 1: Genome → Duplicated Genes ✓
**Status:** Implemented in `pipeline_1/`  
- [x] Extract genome data from Ensembl Plants
- [x] Filter isoforms to get longest protein sequences
- [x] All-vs-all BLAST for similarity detection
- [x] Compute coverage metrics
- [x] Filter protein pairs (id30, qcov50, scov50)
- [x] Build similarity network
- [x] Cluster into gene families using MCL
- **Output:** `output/clusters/protein_families_*.tsv` and `.txt`

### 2. Pipeline Part 2: Duplicated Genes → Ks Values
**Status:** Implemented in `pipeline_2/`  
- [x] Align protein sequences (MUSCLE/MAFFT)
- [x] Back-translate to codon alignment (pal2nal)
- [x] Generate phylip format files
- [x] Calculate Ks values using PAML yn00
- [x] Quality control on Ks values (remove Ks > 5, saturated pairs)
- [x] Output consolidation and validation

### 3. Data Preprocessing & Statistics
- [ ] **Remove mitochondrial and chloroplast genes**
  - Filter gene families to exclude organellar genomes
  - Document chromosome/scaffold filtering criteria
  - Generate clean gene list for downstream analysis
- [x] **Basic statistics:**
  - Total number of genes in genome
  - Number of duplicated genes (in families)
  - Number of singleton genes
  - Distribution of family sizes
  - Coverage statistics per chromosome

---

## Analysis Tasks

### 4. Ks Distribution Analysis
- [x] Generate Ks distribution plots (histograms, density plots)
- [x] Identify peaks corresponding to duplication events
- [x] **Check for WGD (Whole Genome Duplication) events** in Glycine max
  - Literature review: known WGD events in soybean
  - Compare Ks peaks with expected WGD ages
  - Annotate duplication mechanisms (WGD, tandem, proximal, dispersed)
- [x] Age distribution of duplications (convert Ks to millions of years)

### 5. TE (Transposable Elements) Annotation & Analysis
**Critical:** Ensure TE annotation matches genome version (Ensembl Plants Release 41)

#### 5a. TE Data Retrieval & Processing
- [ ] Download TE annotation from Ensembl Plants (GFF format)
  - Verify genome version compatibility
  - Extract TE positions, types, and classifications
- [ ] Parse TE GFF file
  - Classify by TE superfamily (LTR, LINE, SINE, DNA transposons, etc.)
  - Extract TE coordinates (chromosome, start, end, strand)

#### 5b. TE Distribution Analysis
- [ ] **TE numbers and types:**
  - Count TEs by superfamily/family
  - TE length distributions
  - TE age estimation (if LTR divergence data available)
- [ ] **TE genomic distribution:**
  - TE density per chromosome
  - TE distribution relative to genes (proximity analysis)
  - Pericentromeric vs euchromatic regions

#### 5c. TE Coverage Computation
- [ ] **Window-based coverage analysis:**
  - Define bin/window size (e.g., 100kb, 1Mb windows)
  - Calculate TE coverage (bp) per window
  - Calculate TE density (number of TEs) per window
- [ ] **Test TE density in specific regions:**
  - Near duplicated genes vs singleton genes
  - Near TAGs vs non-TAGs
  - Gene-rich vs gene-poor regions
- [ ] Visualization: TE coverage heatmaps/tracks along chromosomes

### 6. TAGs (Tandemly Arrayed Genes) Analysis
**Definition:** Genes in close physical proximity (same chromosome, within X kb)

#### 6a. TAG Identification
- [ ] Define tandem duplication criteria (e.g., ≤10 genes apart, ≤100kb distance)
- [ ] Identify TAG clusters from duplicated gene pairs
- [ ] **Check gene orientation within TAGs:**
  - Same orientation vs opposite
  - Statistical significance of orientation bias
- [ ] Annotate genes as TAG or non-TAG

#### 6b. TAG Age Analysis (Ks-based)
- [ ] Compare Ks distributions: TAGs vs non-TAGs
- [ ] Hypothesis: TAGs are younger than other duplication types
- [ ] Statistical tests (Mann-Whitney U, t-test) for age differences
- [ ] Visualize age distributions (boxplots, violin plots)

#### 6c. TAG Functional Analysis
- [ ] **GO enrichment analysis:**
  - TAGs vs non-TAGs
  - Young TAGs vs old TAGs
- [ ] **Research questions:**
  - Do younger TAGs have specific functions?
  - Is there a link between duplication age and function?
  - Are TAGs enriched for specific biological processes?
- [ ] Tools: clusterProfiler (R), Panther, or other ORA methods

---

## Additional Tasks

### 7. Duplication Type Annotation
- [ ] Classify all duplicated genes by mechanism:
  - WGD (whole genome duplication)
  - Tandem (TAGs)
  - Proximal (nearby but not tandem)
  - Dispersed (different chromosomes or distant)
- [ ] Create annotation file mapping genes to duplication types

### 8. Chromosome-level Analysis
- [ ] Filter scaffolds: focus on whole chromosomes only
- [ ] Create gene mapping file (gene → chromosome → position)
- [ ] Synteny analysis (if comparing to other species)

---

## Extensions & Future Work

- [ ] **Transcriptomics integration:**
  - ArrayExpress dataset E-MTAB-4279 (Glycine max RNA-seq)
  - Expression analysis of duplicated genes
  - TAGs vs non-TAGs expression patterns
  
- [ ] **Comparative genomics:**
  - Ortholog analysis across legume species
  - Conservation of duplication patterns
  
- [ ] **Functional enrichment:**
  - GO enrichment (clusterProfiler)
  - KEGG pathway analysis
  - Domain enrichment

- [ ] **Advanced TAG analysis:**
  - Spatial clustering patterns
  - Orientation statistical tests
  - Expression correlation within TAG clusters

- [ ] **PPI Analysis:**
  - Retrieve PPIs from STRING database
  - Compare PPI networks: TAGs vs non-TAGs
  - Network topology analysis


---

## Data & Documentation

### Key Files
- Genome: Ensembl Plants Release 41
- Protein sequences: `data/protein_info_longest.csv`
- Gene families: `output/clusters/protein_families_*.tsv`
- Ks values: `pipeline_2/` outputs

### Documentation to Create
- [ ] TE analysis methodology document
- [ ] TAG identification criteria and validation
- [ ] Statistical methods for all comparisons
- [ ] Data provenance log (all downloaded files with URLs and versions)

---
