# Duplicated Genes

*currently looking into orientation analysis and circos configuration*

Associated with this part are the following scripts and analysis:
```
├── analysis/
│   └── duplicated_genes/
│       ├── explore_*.R     # -- just exploration no need to run
│       ├── detect_TAGs.R   # -- detect TAGs script
│       ├── load_variables.R
│       ├── plot_*.R        # -- plotting scripts
│       └── duplicated_genes.Rproj
└──scripts/
    ├── compute_TAGs.sh             # -- compute TAGs ratios for spacers 0-10 from detect_TAGs.R
    ├── filter_duplicated_info.sh   # -- filter duplicated genes info out of prot_longest_info.csv
    ├── get_TAGs_ratio.sh           # -- get TAGs counts and ratios for spacers 0-10
    ├── get_TAGs_list.sh            # -- get list of TAGs
    └── get_duplication_ratio.sh    # -- get duplication ratios for different cluster files
```

> [!TIP]
> In case wanted too replicate analyses, recommended to open the .Rproj file in RStudio to have all paths set correctly and run `load_variables.R` first to load all variables used in the analysis scripts

Relevant outputs:
* `output/duplication_classes/TAGs.txt` : contains TAGs list 
* `output/duplication_classes/TAGs/TAGs_1.tsv` : contains TAGs detected with spacer=1 (have info on number of genes in each TAG array, positions, etc)


## Methods to identify duplicated genes

The default thresholds used in our pipeline: 30% identity, 50% coverage => 99% of genes are duplicated genes (89598 isoforms, 56680 genes, 56144 out of them are duplicated genes)

In [^3], the work was performed on rice, with a relatively higher number of initial genes (42534, higher than arabidpsis studies and closer to our number), they used to filtering strategies:  
- 30% identity and 70% coverage => low stringency dataset (reduced to ~18k)
- 50% identity and 90% coverage => high stringency dataset (reduced to ~9k)

> [!TIP]
> maybe would consider the data we have as a low stringency dataset and then perform analysis using a high stringency one (like 70,90) to see how results change

Running pipeline on different thresholds to see how results change, by comparing number of duplicated genes in each file in `output/clusters/` containing gene families detected with different thresholds
```bash
./scripts/get_duplication_ratios.sh
# -- outputs in output/statistics/duplication_ratios.tsv this table
```

| cluster_file | duplicated_genes | total_genes | duplication_ratio |
|--------------|------------------|-------------|-------------------|
| protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv | 56145 | 56679 | 0.9905 |
| protein_families_filtered_blast_results_id30_qcov70_scov70_wcol12_network.tsv | 54271 | 56679 | 0.9575 |
| protein_families_filtered_blast_results_id50_qcov90_scov90_wcol12_network.tsv | 49334 | 56679 | 0.8704 |
| protein_families_filtered_blast_results_id70_qcov90_scov90_wcol12_network.tsv | 44925 | 56679 | 0.7926 |

*consider filtering by # of protein hits, [scripts/top5_per_protein.py](../scripts/top5_per_protein.py) for reference*

**MCScanX** is a widely used tool for detecting and classifying duplicated genes. It provides scripts to classify gene duplicates from blast results and gene position files

> To run MCScanX try to follow these steps or else might give an error:
> * Create a folder called gm/ and put inside it gm.blast (renamed raw blast results with 12 columns) and gm.gff (gff of the shape: chr# start end gene_id)
> * gff taken from prot_longest_info.csv where peptide_id are taken as gene_id (or genomicRanges file taken from R, there s a command to create this file)
> * gff chr# should be of the form gm1, gm2, ..., gm20
> make sure both files are tab seperated and not space seperated
> * run MCScanX: `<path_to_mcscanx>/MCScanX <path_to_gm_dir>/gm`

Results from MCScanX:
```text
$ ./MCScanX gm
Reading BLAST file and pre-processing
Generating BLAST list
0 matches imported (9336543 discarded)
0 pairwise comparisons
0 alignments generated
Pairwise collinear blocks written to gm.collinearity [85.624 seconds elapsed]
Writing multiple syntenic blocks to HTML files
...
Done! [20.653 seconds elapsed]
```
```text
$ ./duplicate_gene_classifier gm
Reading BLAST file and pre-processing
Generating BLAST list
0 matches imported (9336543 discarded)
0 pairwise comparisons
0 alignments generated
Type of dup     Code    Number
Singleton       0       56181
Dispersed       1       0
Proximal        2       0
Tandem  3       0
WGD or segmental        4       0
```

It is also recommended by the tool developers to filter the BLAST results to keep only the top 5 hits per gene in a command like this (tried it but results not too different except for number of discarded matches is much less):

```bash
# tmux new -s mcscanx_blast
makeblastdb -in ../../data/peptides_longest.fa -dbtype prot -out peptides_db
blastp -query ../../data/peptides_longest.fa -db peptides_db   -evalue 1e-10 -max_target_seqs 5 -out gm.blast -outfmt 6
```

## Duplicated genes studies in _Glycine max_

Studies in Glycine max show the organism to have around 70-75% of its genes as duplicated genes [^4][^5], due to its paleopolyploid nature with 2 rounds of whole-genome duplications. Yang et al. (2013) [^5] performed their analysis on genomes taken from Plant Genome Duplication Database and used the MCScan tool to perform syntenic black dedtection, classification of dup genes and downstream analysis.

Most studies we found on duplicated genes in Glycine max were family-focused, i.e. they studied the duplication patterns of specific gene families such as PP2C, WRKY genes.., rather than performing a genome-wide analysis of duplicated genes. Some performed genome-wide analysis using BAC sequences and mapping to identify duplicates.


## Types of duplications

### Tandemly arrayed genes (TAGs)

| Script Path | Description / Plot |
|-------------|--------------------|
| [`analysis/duplicated_genes/detect_TAGs.R`](../analysis/duplicated_genes/detect_TAGs.R) | Detect TAGs with a given spacer number |
| [`scripts/get_TAGs_ratio.sh`](../scripts/get_TAGs_ratio.sh) | Automates detection for spacer numbers 0–10 (results in `output/duplication_classes/TAGs/TAGs_<spacer_number>.tsv`) |
| [`analysis/duplicated_genes/plot_TAGs_distribution.R`](../analysis/duplicated_genes/plot_TAGs_distribution.R) | [![tag_vs_spacer](./assets/TAGs_vs_spacer.png)](../analysis/duplicated_genes/plot_TAGs_distribution.R) |
| [`analysis/duplicated_genes/plot_TAGs_distribution.R`](../analysis/duplicated_genes/plot_TAGs_distribution.R) | [![tag array vs spacer](./assets/TAG_array_vs_spacer.png)](../analysis/duplicated_genes/plot_TAGs_distribution.R) |
| [`analysis/duplicated_genes/plot_TAGs_distribution.R`](../analysis/duplicated_genes/plot_TAGs_distribution.R) | [![tag spacer 1 size distribution](./assets/TAG_sizes_dist.png)](../analysis/duplicated_genes/plot_TAGs_distribution.R) |
| [`analysis/duplicated_genes/plot_chr_distribution.R`](../analysis/duplicated_genes/plot_chr_distribution.R) | [![tag chr distribution](./assets/TAGs_chr_dist.png)](../analysis/duplicated_genes/plot_chr_distribution.R)  |


In this review [^1], Tandemly arrayed genes are defined as genes that are physically close on the chromosome and share high sequence similarity. As we have dealth previously with the definition of "sequence similarity" to be considered duplicated, what's left to define is "physically close". The gene spacer strategy revolves aroung setting a max spacer number, i.e. threshold of intervening genes, to consider two genes as tandemly duplicated. Usually this spacer number ranges between 0 (a perfect TAG cluster with no intervening genes) to 10, with 1, 5 and 10 being common choices. To choose ours, we refer to Shoja & Zhang (2006) [^2] who tried the 11 different spacer numbers from 0 to 11 on 3 different genomes (human, mouse, rat) and observed the increase in the number of TAGs detected with increasing spacer number. 

We tried the same approach on our data, running the detection of TAGs with spacer numbers from 0 to 10 and plotting the results. We observed a similar trend as in[^2], with a rapid increase in the number of TAGs detected from spacer 0 to 1, then a slower one follows. Bsed on the same approach, we will consider spacer number 1 as our threshold to define TAGs in our data, this is from running the script `script/compute_TAGs.sh` which uses the R script `analysis/duplicated_genes/detect_TAGs.R` to detect TAGs with a given spacer number on 0:10 range (can be found in `output/statistics/TAGs_ratios.tsv`)

| spacer | n_TAG_genes | n_TAG_arrays | ratio_TAG_genes |
|--------|-------------|--------------|-----------------|
| 0      | 8326        | 1201         | 0.14689 |
| 1      | 9627        | 1326         | 0.16985 |
| 10     | 11455       | 1441         | 0.20210 |
| 2      | 10150       | 1371         | 0.17907 |
| 3      | 10475       | 1396         | 0.18481 |
| 4      | 10663       | 1408         | 0.18812 |
| 5      | 10826       | 1417         | 0.19100 |
| 6      | 10977       | 1423         | 0.19366 |
| 7      | 11102       | 1427         | 0.19587 |
| 8      | 11244       | 1432         | 0.19838 |
| 9      | 11357       | 1436         | 0.20037 |


[![tag_vs_spacer](./assets/TAGs_vs_spacer.png)](../analysis/duplicated_genes/plot_TAGs_distribution.R) 

So for tags, as spacer 1 is chosen, the file [`output/duplication_classes/TAGs/TAGs_1.tsv`](../output/duplication_classes/TAGs/TAGs_1.tsv) will be used in downstream analysis (seqnames is chr# btw)


### Proximal and dispersed duplications

Proximal can be defined as duplicated genes that are seperated by a few number of genes/kbs on the same chromosome, while dispersed are duplicated genes that are located far away on the same chromosome or on different chromosomes. To consider proximal will first need to define a threshold ideally in kbs distance (unlike tandem where it's more about genes - higher resolution), and then remove from this set the TAGs already detected. The rest that are not TAGs nor proximal will be considered dispersed (would rethink their classification according to WGDs)


### Whole-genome duplications (WGDs)

Based on Ks plot and values, startegy first:
* detect WGD events at peaks
* define Ks ranges for each WGD event
* classify duplicated genes based on their Ks values into WGD-derived or not

## Transposable elements and duplicated genes

It's either about studying transposed duplications or rather focus on richness of TEs around duplicated genes based on coverage and windows or something like that. Check [this article for mme rizzon](https://academic.oup.com/gbe/article/13/5/evab062/6273345#303241768) that split duplicated genes based on their TE richness, also [this "one code to find them all": a perl tool to parse RepeatMasker output files](https://link.springer.com/article/10.1186/1759-8753-5-13#Abs1) if it helps

## Ideas

- [ ] try circos to plot duplicated genes links on chromosomes (might be a lot of links => consider largest family or only TAG links for example, or show distribution in heatmap form)

[^1]: Lallemand, T., Leduc, M., Landès, C., Rizzon, C., & Lerat, E. (2020). An overview of duplicated gene detection methods: why the duplication mechanism has to be accounted for in their choice. Genes, 11(9), 1046.  
[^2]: Shoja, V., & Zhang, L. (2006). A roadmap of tandemly arrayed genes in the genomes of human, mouse, and rat. Molecular biology and evolution, 23(11), 2134-2141.
[^3]: Lallemand, T., Leduc, M., Landès, C., Rizzon, C., & Lerat, E. (2020). An overview of duplicated gene detection methods: why the duplication mechanism has to be accounted for in their choice. Genes, 11(9), 1046.
[^4]: Kim, K. D., El Baidouri, M., Abernathy, B., Iwata-Otsubo, A., Chavarro, C., Gonzales, M., ... & Jackson, S. A. (2015). A comparative epigenomic analysis of polyploidy-derived genes in soybean and common bean. Plant Physiology, 168(4), 1433-1447.
[^5]: Yang, Y., Wang, J., & Di, J. (2013). Comparative inference of duplicated genes produced by polyploidization in soybean genome. International journal of genomics, 2013(1), 275616.