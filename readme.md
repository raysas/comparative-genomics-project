# comparative-genomics-project

Hey, best place to put your code is in analysis, if you can make a folder called `family_sizes` or smtg unqiue to your analysis. Kindly find the output files in the output folder

At this stage we have 4 main new outputs, where almost each will comprise 2 outputs:

- "low" or "id30_cov50_evalue1e-10" named: low stringency dataset at filtration of 30% identity, 50% coverage and evalue 1e-10
- "high" or "id50_cov70_evalue1e-10" named: high stringency dataset at filtration of 50% identity, 70% coverage and evalue 1e-10

*might be worth running analysis on each of them*

This is a description of each, i assume you might need gene_lists/ the most, to get full list of genes (dups and non dups, u can use output/info/protein_longest_info.csv to get full info on all genes including non duplicated one or create non tag duplicated genes from there as well)*

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