# TODO - Comparative Genomics Project

## Major Themes

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
Family sizes analysis:

- [ ] functional annotation of gene families based on sizes (seperate to small, and large families + singletons maybe)
- [x] family size distribution plots
- [ ] largest family analysis (functional, ppi, identify what family, distribution on genome, etc)

> [!NOTE]
> will have both high stringency (H) and low stringency (L) results which will lead to different families in teh results, the analysis most probably would be related to the high stringency one

General tasks:

- [x] pipeline 1: from raw data to duplicated genes
- [x] pipeline 2: from duplicated genes to Ks computation
- [ ] annotate duplicated genes with different types (TAGs, WGD...)
- [ ] retrieval of PPIs from STRING database from a list of genes
- [ ] General analysis of duplicated genes families
- [ ] General analysis of Ks distributions and age of duplications
=======
- [ ] pipeline 1: from raw data to duplicated genes
- [ ] pipeline 2: from duplicated genes to Ks computation
- [ ] duplicated genes classification (TAGs, WGD...)
- [ ] comparative analysis of singletons and duplicated genes families
- [ ] analysis of Ks distributions and age of duplications
>>>>>>> master
- [ ] TEs annotation and analysis
- [ ] Computation of TE coverage


### Extensions

- [ ] Transcriptomics data integration into the project: check [this dataset](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-4279?query=glycine%20max) from arrayExpress for RNA-seq data of Glycine max - accession:E-MTAB-4279
- [ ] Comparative orthologs analysis: based on fully reproducible one-liner pipeline integration
- [ ] Functional analysis of duplicated genes: at least GO enrichment (could be on panther or using R packages like clusterProfiler, it should be following the same ORA method based on statistical tests)
- [ ] Localization and orientation of TAGs (also can do some functional/ppi types of analysis to compare)
<<<<<<< HEAD
=======
## Duplicated Genes

- [ ] functional analysis related to Ks and Ka and age
- [ ] functional analysis in relation to duplication type
- [ ] functional analysis related to family sizes
- [ ] WGD genes retrieval and analysis
- [ ] MCScanX and circos
>>>>>>> general_analysis
=======
- [ ] TE annotation
- [ ] Comaprison between databases
- [ ] Coverage computation
- [ ] TE landscape plots (?)
- [ ] Integration with duplicated genes and Ks plots (??)
<<<<<<< HEAD
>>>>>>> TE
=======


Joelle


- [ ] Proportion of TEs in different genomic regions (e.g., upstream, downstream, intronic, intergenic)
- [ ] Correlation analysis between TE density and gene density
- [ ] Age distribution of TEs based on divergence from consensus sequences
- [ ] Identification of recent TE insertions
- [ ] Analysis of TE impact on gene expression (if RNA-seq data is available)
- [ ] Comparative analysis of TE content across different species or strains
- [ ] Visualization of TE landscapes using tools like RepeatMasker or custom scripts
- [ ] Functional annotation of genes disrupted by TE insertions


More:
- [] DIstribution plots 
- [] percentages of TEs in different regions (upstream, downstream, intronic, intergenic)
- [] TE density across the genome (in each chromosome)
  



### **Project Scope: Integrating Gene Duplication and TE Dynamics in *Glycine max***

This project aims to reconstruct the evolutionary history of the *Glycine max* genome by integrating two parallel analyses: **Whole Genome Duplication (WGD)** and **Transposable Element (TE) proliferation**. We will analyze the relationship between these features across three dimensions: Time, Space, and Genomic Architecture.

---

### **Phase 1: Independent Track Analysis**

1.  **The Gene Track (WGD Events):**
    * Identify paralogous gene pairs (duplicated genes).
    * Calculate Synonymous Substitution rates ($K_s$) for these pairs.
    * **Goal:** Identify the $K_s$ peak corresponding to the *Glycine* WGD event (~13 MYA).

2.  **The TE Track (Transposition Bursts):**
    * Map TEs using the Atlas reference library.
    * Calculate divergence ($K_s$) of TEs from their consensus sequences.
    * **Goal:** Identify peaks in TE accumulation history (bursts).

---

### **Phase 2: Integrative Analysis (The Synthesis)**

#### **1. Temporal Integration: The "Genomic Shock" Hypothesis**
* **Method:** Overlay the Gene $K_s$ distribution and the TE $K_s$ distribution on a single timeline.
* **The Question:** Did the Whole Genome Duplication trigger a proliferation of TEs?
* **Hypothesis:** **Genomic Shock.** The stress of polyploidization may have disrupted epigenetic silencing mechanisms, leading to a massive TE burst immediately following the WGD event (overlapping or sequential $K_s$ peaks).

#### **2. Spatial Integration: The "Real Estate" Analysis**
* **Method:** Calculate physical distances between TEs and specific gene categories (WGD-derived duplicates vs. Singletons).
* **The Question:** How does TE proximity influence gene fate?
* **Hypothesis A: Gene Protection vs. Pseudogenization.**
    * *Preserved Duplicates:* Likely reside in TE-poor regions (purifying selection protects them from insertion).
    * *Singletons:* We hypothesize that the "lost" copy of a gene pair may have been disrupted by TE insertion, explaining its deletion or conversion to a pseudogene.
* **Hypothesis B: Neofunctionalization via Promoters.**
    * Analyze TEs inserted in **promoter regions** vs. **introns**.
    * TEs in promoters can alter expression patterns, potentially driving the divergence of one gene copy to acquire a new function (Neofunctionalization) while the other retains the ancestral role.

#### **3. Architectural Integration: Sub-genome Fractionation**
* **Method:** Compare TE density across the two putative sub-genomes resulting from the paleopolyploidy.
* **The Question:** Is TE accumulation biased toward one set of chromosomes?
* **Hypothesis:** **Biased Fractionation.** In paleopolyploids, one sub-genome usually remains dominant while the other degrades. We hypothesize that the degrading (fractionated) sub-genome will show a significantly higher density of TEs, serving as a marker for genomic decay.
>>>>>>> TE
=======
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
