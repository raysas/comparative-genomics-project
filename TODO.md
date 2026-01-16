# TODO

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
- [ ] TEs annotation and analysis
- [ ] Computation of TE coverage

Extensions:

- [ ] Transcriptomics data integration into the project: check [this dataset](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-4279?query=glycine%20max) from arrayExpress for RNA-seq data of Glycine max - accession:E-MTAB-4279
- [ ] Comparative orthologs analysis: based on fully reproducible one-liner pipeline integration
- [ ] Functional analysis of duplicated genes: at least GO enrichment (could be on panther or using R packages like clusterProfiler, it should be following the same ORA method based on statistical tests)
- [ ] Localization and orientation of TAGs (also can do some functional/ppi types of analysis to compare)
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
