# Duplicated genes

# On Glycine max

Studies in Glycine max show the organism to have around 70-75% of its genes as duplicated genes [^4][^5], due to its paleopolyploid nature with 2 rounds of whole-genome duplications. Yang et al. (2013) [^5] performed their analysis on genomes taken from Plant Genome Duplication Database and used the MCScan tool to perform syntenic black dedtection, classification of dup genes and downstream analysis.

Most studies we found on duplicated genes in Glycine max were family-focused, i.e. they studied the duplication patterns of specific gene families such as PP2C, WRKY genes.., rather than performing a genome-wide analysis of duplicated genes. Some performed genome-wide analysis using BAC sequences and mapping to identify duplicates.

## On pipeline 1

Starting from blast results, we identify duplicated genes in the genome.

Done some diagnostics plots on the blast results to see how identity and coverage are distributed (script:[`scripts/blast_diagnostics.py`](../scripts/blast_diagnostics.py))


| id distribution | coverage distribution | id vs coverage |
|-----------------|-----------------------|----------------|
| ![id distribution](./assets/blast_viz_pident_hist.png) | ![coverage distribution](./assets/blast_viz_coverage_hist.png) | ![id vs coverage](./assets/blast_viz_pident_vs_cov.png) |

The default thresholds to rub blastp used in our pipeline: 30% identity, 50% coverage => 99% of genes are duplicated genes (89598 isoforms filtered to 56680 genes, 56144 out of them are duplicated genes)

In [^3], the work was performed on rice, with a relatively higher number of initial genes (42534, higher than arabidpsis studies and closer to our number), they used to filtering strategies:  
- 30% identity and 70% coverage => low stringency dataset (reduced to ~18k)
- 50% identity and 90% coverage => high stringency dataset (reduced to ~9k)


## On datasets

At the start, on an exploratory phase, we tried performing the full analysis on a low stringency dataset (30% id, 50% cov) to check initial resultsm knowing that we can use a more stringent dataset later on based on filtration thresholds. The dataset had a 99% duplication rate which led to a pretty exhaustive search space to work with in pipeline 2 in particular, requiring really high computational power to run (yazid part here).
From a large search space of the righ filtration thresholds, we gor preliminary statistics on duplicated genes as well as some analysis results that made us judge better what to go for, for the rest of the analysis.  
Running pipeline on different thresholds to see how results change, by comparing number of duplicated genes in each file in `output/clusters/` containing gene families detected with different thresholds
```bash
./scripts/dups_get_ratio.sh
# -- outputs in output/statistics/duplication_ratios.tsv this table
```

| id  | qcov | scov | evalue | duplicated_genes | total_genes | duplication_percent | largest_family_size |
|-----|------|------|--------|------------------|-------------|----------------------|----------------------|
| 30  | 50   | 50   | 1e-10 | 52139 | 56679 | 91.98% | 1326 |
| 30  | 50   | 50   | 1e-10 | 52031 | 56679 | 91.79% | 1323 |
| 30  | 50   | 50   | 1e-5  | 52516 | 56679 | 92.65% | 1311 |
| 30  | 50   | 50   | none  | 56145 | 56679 | 99.05% | 1314 |
| 30  | 60   | 60   | 1e-10 | 51870 | 56679 | 91.51% | 1324 |
| 30  | 60   | 60   | 1e-10 | 51768 | 56679 | 91.33% | 1321 |
| 30  | 60   | 60   | 1e-5  | 52227 | 56679 | 92.14% | 1309 |
| 30  | 70   | 70   | 1e-10 | 51196 | 56679 | 90.32% | 1320 |
| 30  | 70   | 70   | 1e-10 | 51114 | 56679 | 90.18% | 1302 |
| 30  | 70   | 70   | 1e-5  | 51526 | 56679 | 90.90% | 1304 |
| 30  | 70   | 70   | none  | 54271 | 56679 | 95.75% | 1307 |
| 40  | 60   | 60   | 1e-10 | 51484 | 56679 | 90.83% | 828 |
| 40  | 60   | 60   | 1e-5  | 51814 | 56679 | 91.41% | 842 |
| 40  | 70   | 70   | 1e-10 | 50801 | 56679 | 89.62% | 825 |
| 40  | 70   | 70   | 1e-5  | 51101 | 56679 | 90.15% | 840 |
| 40  | 80   | 80   | 1e-10 | 49824 | 56679 | 87.90% | 821 |
| 40  | 80   | 80   | 1e-5  | 50085 | 56679 | 88.36% | 835 |
| 50  | 70   | 70   | 1e-10 | 50031 | 56679 | 88.27% | 350 |
| 50  | 70   | 70   | 1e-10 | 49956 | 56679 | 88.13% | 349 |
| 50  | 70   | 70   | 1e-5  | 50280 | 56679 | 88.71% | 357 |
| 50  | 80   | 80   | 1e-10 | 49073 | 56679 | 86.58% | 349 |
| 50  | 80   | 80   | 1e-5  | 49290 | 56679 | 86.96% | 355 |
| 50  | 90   | 90   | 1e-10 | 48080 | 56679 | 84.82% | 348 |
| 50  | 90   | 90   | 1e-5  | 48261 | 56679 | 85.14% | 354 |
| 50  | 90   | 90   | none  | 49334 | 56679 | 87.04% | 360 |
| 60  | 70   | 70   | 1e-10 | 48662 | 56679 | 85.85% | 211 |
| 60  | 70   | 70   | 1e-5  | 48863 | 56679 | 86.21% | 214 |
| 60  | 80   | 80   | 1e-10 | 47723 | 56679 | 84.19% | 209 |
| 60  | 80   | 80   | 1e-5  | 47897 | 56679 | 84.50% | 212 |
| 70  | 90   | 90   | 1e-10 | 44341 | 56679 | 78.23% | 177 |
| 70  | 90   | 90   | none  | 44925 | 56679 | 79.26% | 177 |

We chose 2 datasets to work with for the rest of the analysis:
- Low stringency: 30% id, 50% cov, 1e-10 evalue => 91.98% duplication rate
- High stringency: 50% id, 70% cov, 1e-10 evalue => 88.27% duplication rate

Main results for each can be found in files suffixed with `_low` and `_high` respectively in subdirectoried of the `output/` folder.

## On gene families

*gene family size distribution*

### largest family

## On TAGs

### gene spacer based indentification

The first approach used to identify tandemly arrayed genes (TAGs) was based on several works like [^1], also mentioned in this review [^2].
In the review [^1], Tandemly arrayed genes are defined as genes that are physically close on the chromosome and share high sequence similarity. As we have dealth previously with the definition of "sequence similarity" to be considered duplicated, what's left to define is "physically close". The gene spacer strategy revolves aroung setting a max spacer number, i.e. threshold of intervening genes, to consider two genes as tandemly duplicated. Usually this spacer number ranges between 0 (a perfect TAG cluster with no intervening genes) to 10, with 1, 5 and 10 being common choices. To choose ours, we refer to Shoja & Zhang (2006) [^2] who tried the 11 different spacer numbers from 0 to 10 on 3 different genomes (human, mouse, rat) and observed the increase in the number of TAGs detected with increasing spacer number. 

We tried the same approach on our data, running the detection of TAGs with spacer numbers from 0 to 10 and plotting the results. We observed a similar trend as in[^2], with a rapid increase in the number of TAGs detected from spacer 0 to 1, then a slower one follows. Bsed on the same approach, we will consider spacer number 1 as our threshold to define TAGs in our data, this is from running the script `script/compute_TAGs.sh` which uses the R script `scripts/detect_TAGs.R` to detect TAGs with a given spacer number on 0:10 range.

```bash
Rscript scripts/TAGs_detect.R --help
./scripts/TAGs_compute.sh -h
# Usage: ./scripts/TAGs_compute.sh [low|high]
```

* low stringency dataset ([`output/statistics/TAGs_spacers_ratios_low.tsv`](../output/statistics/TAGs_spacers_ratios_low.tsv) )

| spacer | n_TAG_genes | n_TAG_arrays | TAG_genes_percent |
|--------|-------------|--------------|--------------------|
| 0      | 8563        | 1263         | 15.107% |
| 1      | 9696        | 1381         | 17.106% |
| 2      | 10164       | 1425         | 17.932% |
| 3      | 10442       | 1447         | 18.423% |
| 4      | 10607       | 1458         | 18.714% |
| 5      | 10751       | 1468         | 18.968% |
| 6      | 10894       | 1474         | 19.220% |
| 7      | 11010       | 1480         | 19.425% |
| 8      | 11151       | 1487         | 19.673% |
| 9      | 11267       | 1495         | 19.878% |
| 10     | 11370       | 1499         | 20.060% |

* high stringency dataset ([`output/statistics/TAGs_spacers_ratios_high.tsv`](../output/statistics/TAGs_spacers_ratios_high.tsv) )

| spacer | n_TAG_genes | n_TAG_arrays | TAG_genes_percent |
|--------|-------------|--------------|--------------------|
| 0      | 7524        | 1597         | 13.274% |
| 1      | 8545        | 1756         | 15.076% |
| 2      | 8926        | 1821         | 15.748% |
| 3      | 9122        | 1852         | 16.094% |
| 4      | 9239        | 1863         | 16.300% |
| 5      | 9321        | 1875         | 16.445% |
| 6      | 9393        | 1886         | 16.572% |
| 7      | 9447        | 1892         | 16.667% |
| 8      | 9501        | 1898         | 16.762% |
| 9      | 9551        | 1908         | 16.851% |
| 10     | 9605        | 1917         | 16.946% |


Both datasets showed the same trends with %s being more or less similar:


| low | high |
|-----|------|
| ![](./assets/TAGs_spacers_low.png) | ![](./assets/TAGs_spacers_high.png) |
| ![](./assets/TAGs_array_spacers_high.png) | ![](./assets/TAGs_array_spacers_low.png) |

So for TAGs, gene spacer=1 is chosen as threshold for both datasets. Looking closely at them we see:

| low | high |
|-----|------|
| ![](./assets/TAGs_size_low.png) | ![](./assets/TAGs_size_high.png) |
| ![](./assets/TAGs_chr_low.png) | ![](./assets/TAGs_chr_high.png) | 

same distribution in general, if u dont want to show it, maybe eya or nhi can show their part (if relevant) as preperation to justify that functional analysis for the low and high will overall give similar results

### yazid mcscan and kbs part

## Orientation of TAG gene pairs

*performed on spacer results*

TAGs can be classified into three categories based on the relative orientation of the genes in the pair: tandem (→ → or ← ←), convergent (→ ←), and divergent (← →). The orientation of TAG gene pairs can provide insights into their evolutionary history and functional relationships.  
As TAGs originate from duplicated genes that are positioned next to each other on a chromosome, and unequal crossing-over (UCO) is a recombination-based process in which misaligned homologous chromosomes exchange genetic material, generating these tandem duplicates typically in direct orientation.
According to Rizzon et al (2006) [^6], tandemly arrayed genes (TAGs) are predominantly arranged in direct orientation, reflecting the influence of unequal crossing-over (UCO) as a major mechanism generating TAGs. It showed that around 80% of TAGs in rice and 88% in Arabidopsis are in direct orientation, suggesting that UCO is likely the primary driver of tandem gene duplication, although other evolutionary processes such as selection, recombination, gene gain, and loss may also play a role.  
In fact, Shoja & Zhang (2006) [^2] also observed that in human, mouse, and rat genomes, the majority of TAGs are arranged in direct orientation, with percentages being 68%, 76%, and 72%, respectively. And they performed a chi-squared test to confirm the significance of these percentages compared to what would be expected by chance.

Here, we want to analyze the orientation of out TAG gene pairs to see if they follow the same trend as in [^6]. 

```bash
Rscript scripts/TAGs_pairs_orientation.R --stringency low
Rscript scripts/TAGs_pairs_orientation.R --stringency high
```

* low stringency dataset
    ```
        orientation count percentage
    <chr>       <int>      <dbl>
    1 tandem      14425      82.5 
    2 convergent   1637       9.36
    3 divergent    1423       8.14

            Chi-squared test for given probabilities

    data:  observed
    X-squared = 7392.3, df = 2, p-value < 2.2e-16
    ```
* high stringency dataset
    ```
    orientation count percentage
    <chr>       <int>      <dbl>
    1 tandem      12685      83.5 
    2 convergent   1320       8.69
    3 divergent    1192       7.84

            Chi-squared test for given probabilities

    data:  observed
    X-squared = 6812, df = 2, p-value < 2.2e-16
    ```

*table and results copied from Rscript output*  
The extremely low p-value indicates that the observed distribution of TAG gene pair orientations is significantly different from what would be expected by chance, suggesting a strong bias towards tandem orientation among TAGs in this dataset.


[^1]: Lallemand, T., Leduc, M., Landès, C., Rizzon, C., & Lerat, E. (2020). An overview of duplicated gene detection methods: why the duplication mechanism has to be accounted for in their choice. Genes, 11(9), 1046.  
[^2]: Shoja, V., & Zhang, L. (2006). A roadmap of tandemly arrayed genes in the genomes of human, mouse, and rat. Molecular biology and evolution, 23(11), 2134-2141.
[^3]: Lallemand, T., Leduc, M., Landès, C., Rizzon, C., & Lerat, E. (2020). An overview of duplicated gene detection methods: why the duplication mechanism has to be accounted for in their choice. Genes, 11(9), 1046.
[^4]: Kim, K. D., El Baidouri, M., Abernathy, B., Iwata-Otsubo, A., Chavarro, C., Gonzales, M., ... & Jackson, S. A. (2015). A comparative epigenomic analysis of polyploidy-derived genes in soybean and common bean. Plant Physiology, 168(4), 1433-1447.
[^5]: Yang, Y., Wang, J., & Di, J. (2013). Comparative inference of duplicated genes produced by polyploidization in soybean genome. International journal of genomics, 2013(1), 275616.
[^6]: Rizzon, C., Ponger, L., & Gaut, B. S. (2006). Striking similarities in the genomic distribution of tandemly arrayed genes in Arabidopsis and rice. PLoS computational biology, 2(9), e115.

<!-- "Genes without hits that met a threshold of e-value 10e-10 were deemed singletons" -->