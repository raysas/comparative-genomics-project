# comparative-genomics-project

Things done but still documenting and polishing:
- [x] tag spacer analysis + plots for low and high (shows same pattern and ~%s)
- [x] orientation analysis 
- [ ] in the process: looking at biggest array of TAGs - function and viz
- [x] largest family functional analysis (seem like a dna helicase family)
- [ ] in the making: biggest cluster in network viz
- [ ] refixed codes on new dir structure
- [ ] ..........

*commiting to checkout and move some things, will fix morning*

At this stage we have 4 main new outputs, where almost each will comprise 2 outputs:

- "low" or "id30_cov50_evalue1e-10" named: low stringency dataset at filtration of 30% identity, 50% coverage and evalue 1e-10
- "high" or "id50_cov70_evalue1e-10" named: high stringency dataset at filtration of 50% identity, 70% coverage and evalue 1e-10

*might be worth running analysis on each of them*


1. [`output/gene_lists`](output/gene_lists): contains lists of genes classified in groups, it has them in txt files where each line is a peptide id, the groups are:
```
output/gene_lists/
├── TAGs
│   └── spacer_based
│       ├── TAGs_high.txt
│       └── TAGs_low.txt
├── archive
│   ├── H_largest_family.txt
│   ├── TAGs.txt
│   ├── largest_family.txt
│   └── singleton_genes.txt
├── largest_family
│   ├── largest_family_high.txt
│   ├── largest_family_id30_cov50_evalue1e-10.txt
│   ├── largest_family_id50_cov70_evalue1e-10.txt
│   └── largest_family_low.txt
├── mapped_peptide_gene_ids
│   ├── largest_family_id30_cov50_evalue1e-10_mapped.txt
│   ├── largest_family_id50_cov70_evalue1e-10_mapped.txt
│   ├── singletons_high_mapped.txt
│   └── singletons_low_mapped.txt
└── singletons
    ├── singletons_high.txt
    └── singletons_low.txt
```


2. [`output/statistics`](output/statistics): contains statistics files, including the ratios of TAG genes and arrays per spacer size, in `TAGs_ratios_*.tsv` files and summary statistics in `duplication_statistics_*.tsv` files on yield of runs on different thresholds

```
output/statistics/
├── TAGs_spacers_ratios_high.tsv
├── TAGs_spacers_ratios_low.tsv
└── duplication_ratios.tsv
```

*will put tables in md in final version*

1. [`output/duplication_classes/`](output/duplication_classes/): contains the main outputs of the identification of duplication types, including TAGs so far through the gene spacer based method, in `TAGs/low` and `TAGs/high` folder, with files named `TAGs_*.tsv` that correspond to different gene spcacer thresholds

```
output/duplication_classes/
├── TAGs
│   ├── high
│   │   ├── TAG_gene_pairs.tsv
│   │   ├── TAGs_0.tsv
│   │   ├── TAGs_1.tsv
│   │   ├── TAGs_10.tsv
│   │   ├── TAGs_2.tsv
│   │   ├── TAGs_3.tsv
│   │   ├── TAGs_4.tsv
│   │   ├── TAGs_5.tsv
│   │   ├── TAGs_6.tsv
│   │   ├── TAGs_7.tsv
│   │   ├── TAGs_8.tsv
│   │   └── TAGs_9.tsv
│   └── low
│       ├── TAG_gene_pairs.tsv
│       ├── TAGs_0.tsv
│       ├── TAGs_1.tsv
│       ├── TAGs_10.tsv
│       ├── TAGs_2.tsv
│       ├── TAGs_3.tsv
│       ├── TAGs_4.tsv
│       ├── TAGs_5.tsv
│       ├── TAGs_6.tsv
│       ├── TAGs_7.tsv
│       ├── TAGs_8.tsv
│       └── TAGs_9.tsv
└── singletons
    ├── singletons_high.csv
    └── singletons_low.csv
```

4. [`output/info`](output/info): contains info files used in the process, including the duplicated genes info files at different stringency levels, subsetted from `protein_longest_info.csv` which is also copied here for reference

```
output/info/
├── duplicated_genes_info_high.csv
├── duplicated_genes_info_id30_qcov50_scov50_evalue1e-10_wcol12.csv
├── duplicated_genes_info_id50_qcov70_scov70_evalue1e-10_wcol12.csv
├── duplicated_genes_info_low.csv
└── protein_info_longest.csv
```

note: duplicated_genes_info_id30_qcov50_scov50_evalue1e-10_wcol12.csv=duplicated_genes_info_low.csv and duplicated_genes_info_id50_qcov70_scov70_evalue1e-10_wcol12.csv=duplicated_genes_info_high.csv


## code

will update tmrw