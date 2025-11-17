# TODO

- [ ] pipeline 1: from raw data to duplicated genes
- [ ] pipeline 2: from duplicated genes to Ks computation
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