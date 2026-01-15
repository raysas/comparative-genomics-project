# Duplication Events

We executed pipeline 2 to produce the pairwise synonymous substitution estimates (Ks) from BLAST-derived homologs (BLAST prefilters: ≥30% identity, ≥50% coverage). These Ks results were then filtered to retain reliable values for duplication-age analysis. Primary filters removed non-numeric or missing Ks, non-positive values (Ks ≤ 0) and saturated/high values (commonly Ks > 2.0) that confound peak detection. Larger Ks values (e.g., >0.75) are associated with increasingly large error due to substitutional saturation and the limitations of synonymous site estimation (Li, 1997). To minimize this error while retaining a reasonably sized data set, we adopted a Ks cutoff of 2.0, following the approach of Blanc and Wolfe (2004). This threshold ensures that the retained gene pairs represent reliable duplication events, while reducing the confounding effects of highly diverged or saturated pairs. Additional considered filters included minimum aligned length, minimum percent identity on the alignment used for Ks, and removal of problematic back-translations (frameshifts/stop codons). At each step we retained both the filtered and the full archived datasets.

Filtered Ks values are then visualized as histograms and smoothed density curves to reveal duplication-age peaks. When converting Ks to approximate ages, we use the molecular clock rate λ = 6.1×10⁻⁹ substitutions/site/year to compute age with the formula: Age (MYA) = Ks / (2 × λ) / 1e6. This lambda rate is reported for soybean and used in multiple publications in the literature (Duan et al., 2023; Miura et al., 2008).

![](plots/ks_plot.png) 

**Figure 1.** Ks distribution of all BLAST-derived homologous gene pairs in soybean, showing major duplication-age peaks.

We note two peaks, one around 0.1 and another around 0.5. This is in accordance with previous publications on soybean. Yang et al. (2013) reported peaks at 0.15 and 0.42, while Rulin et al. (2013) reported them at 0.1 and 0.53. However, the corresponding ages of these two duplication events of glycine max have been reported to be 10 and 56.5 Mya (Kim et al., 2015), or 13 and 59 Mya (Schmutz et al., 2010). Accorindgly, the age of the second peak does not match the ks peak value using the previously highlighted lambda value for conversion. According to Schmutz et al. (2010), this discrepancy is explainable by the fact that the older Ks peak (reported at Ks≈0.59) was assigned to an early-legume WGD by anchoring it to fossil evidence that dates the origin of the papilionoid legumes to ~58–60 Mya, which was then used to reverse compute an effective synonymous substitution rate of ~5.17 × 10⁻⁹ substitutions per site per year. The authors stated: "If the older duplication is assumed to have occurred around 58 Myr ago, then the calculated rate of silent mutations extending back to the duplication would be 5.17 × 10-3".

The more recent Glycine-specific WGD was dated independently using a higher lineage-specific substitution rate, resulting in an estimated age of ~13 Mya. The 6.1×10⁻⁹ lambda could be used for the recent event (glycin max-specific), but another rate has to be introduced to account for the older early-legume event. In our case, we deduced a lambda ~4.1 × 10⁻⁹. We explain the discrepancy from that reported by Schmutz et al. (2010) in the fact that they limited their analysis to families of sizes 2 to 6. 

Rulin et al. also reported a third highly diffuse peak ~1.5 (Schmutz et al., 2010). This is also in accordance with our findings as could be seen in the plot. This likely corresponds to an ancient eudicot whole-genome triplication, which occurred near the origin of core eudicots (~115–130 Mya).

## Anchoring

To further clarify our peaks more distinctively, we used MCScanX to identify collinear anchor pairs (gene pairs retained in conserved blocks across the soybean genome). MCScanX detects these anchors by scanning for regions where gene order and orientation are preserved, indicating large-scale duplication events such as WGDs. By plotting the age distribution of these anchored pairs, we observe much more distinct peaks aligning with the previous reports in the literature, as shown in Fig. This approach ensures that the Ks peaks are specifically associated with large-scale, collinear duplications, providing a robust link between the observed molecular signatures and historical genome duplication events.

![](plots/ks_anchors.png)

**Figure 2.** Ks distribution of collinear anchor pairs identified by MCScanX, highlighting distinct peaks associated with large-scale duplication events.

## Filtering

Using our computed Ks values for all 2.4M pairs (at identity >30% and coverage 50%), we experimented with stricter filtering criteria to examine tge effect on the Ks distribution plot. As highlighted in Fig, we note that only ~12.5% of pairs are retained for Ks range 0-2, with more stringent filtering primarily affecting the huge jump in the 1.5-2 region.

![](plots/ks_comparison_grid.png)

**Figure 3.** Comparison of Ks distributions under different filtering criteria, illustrating the effect of stricter thresholds on the retention of gene pairs.

# Automation

Beyond implementing pipeline 1 and 2 for Ks estimation, we fully automated the entire workflow to enable robust, reproducible, and scalable comparative genomics analyses. Each pipeline step—from data extraction and BLAST searches to alignment, Ks calculation, and downstream filtering—was scripted for end-to-end execution, minimizing manual intervention and reducing error. Batch processing, parallelization, and checkpointing were integrated to efficiently handle millions of gene pairs, optimize runtime, and ensure resilience to interruptions. All parameters (e.g., identity, coverage, Ks thresholds) are configurable, allowing rapid adaptation to new species or datasets and facilitating comparative analysis across genomes. The modular design supports easy replication and extension, with outputs and intermediate results systematically logged for transparency. This automation not only streamlines analysis but also enables efficient exploration of different filtering criteria and thresholds, ensuring that the workflow remains flexible and scalable for future studies. We validated our pipeline on additional species, demonstrating rapid and reliable replication of all analysis steps and outputs for new genomic datasets.

# Duplication Types

## TAGs

We used multiple approaches to obtain TAGs. The venn diagram in Fig highlights the overlap of predictions.
![](plots/tags_overlap_venn.png)

**Figure 4.** Venn diagram showing the overlap of tandem adjacent gene (TAG) predictions from distance-based, genecount-based, and MCScanX approaches.

The discrepancies observed between TAGs identified by the distance-based, genecount-based, and MCScanX approaches is due to their different definitions of tandem adjacency. The distance-based approach classifies TAGs as gene pairs located within a specified physical distance (≤100,000 base pairs by default) on the chromosome, regardless of the number of intervening genes. In contrast, the genecount-based method considers pairs as TAGs if they are separated by no more than 10 genes (by default), independent of the actual base pair distance. MCScanX predicts tandem adjacent genes (TAGs) using a collinearity-based approach. First, MCScanX identifies homologous gene pairs through sequence similarity (typically from BLAST results). It then scans chromosomes for collinear blocks—regions where homologous genes are arranged in the same order and orientation. Within these blocks, MCScanX applies a gap threshold (default: ≤10 intervening genes) to define whether genes are sufficiently close to be considered part of the same block. For TAG prediction specifically, MCScanX looks for homologous gene pairs that are directly adjacent (no intervening genes) on the same chromosome and belong to the same collinear block. Only pairs meeting these strict criteria—direct adjacency, shared block membership, and sequence homology—are classified as TAGs by MCScanX. This approach is more conservative than simple distance or genecount methods, as it requires both physical proximity and evidence of conserved genomic context.

The overlap in the Venn diagram represents TAGs consistently detected by all three methods, while the unique regions highlight pairs identified only by one or two approaches due to these differing criteria. This comparison underscores the importance of method selection and parameter choice in tandem duplication analysis, as each approach captures distinct aspects of genomic organization.

![](plots/ks_TAG_vs_nonTAG_density_2.png)

**Figure 5.** Ks density distributions for TAG versus non-TAG gene pairs, demonstrating the age bias of tandem duplications.

We examined the age distribution of tandem adjacent gene (TAG) pairs, as determined by their Ks values. The distribution reveals a concentration of younger TAGs, reflecting ongoing tandem duplication activity, with a tail of older events. This pattern supports the biological expectation that tandem duplications are a continuous process, contributing to gene family expansion and functional diversification. 

## Structural Analysis

Our structural analysis of duplicated gene pairs in soybean reveals a significant enrichment near chromosome ends, inaccordance with what was reported by Yang et al. (2013). As shown in the bar plot and histogram figures, a substantial proportion of both tandem (TAG) and whole-genome duplication (WGD) pairs are located within 4 Mb of the chromosome termini. Quantitatively, we observe that approximately 40% of duplicated gene pairs are found in these subtelomeric regions. This spatial bias is evident in both the overall genome-wide histogram and the per-chromosome analysis, where the frequency of gene pairs declines with increasing distance from the ends and peaks within the first few megabases.

The per-chromosome histograms further illustrate that this enrichment is consistent across individual chromosomes, with the majority of TAG and WGD pairs clustering near the telomeric regions. The centromeric regions, highlighted in purple, show a marked depletion of duplicated pairs, supporting the notion that chromosomal ends are hotspots for gene duplication retention. These results reinforce the conclusion that the chromosomal environment, particularly proximity to telomeres, plays a significant role in the retention and distribution of duplicated genes in soybean, as previously described by Yang et al. (2013).

![](plots/distance_to_end_chr1_histogram.png)
![](plots/distance_to_end_histogram.png)
![](plots/distance_to_end_percentage.png)
![](plots/zone_counts.png)

**Figure 6.** Distribution of gene pair counts by chromosomal zone (quintiles), for TAG and WGD pairs.

**Figure 7.** Percentage of TAG and WGD gene pairs located within 4 Mb of chromosome ends, indicating enrichment near telomeres.

**Figure 8.** Histogram of distances to chromosome ends for TAG and WGD pairs, with telomere threshold (4 Mb) marked in red and expected centromere region shaded in purple.

**Figure 9.** Per-chromosome histogram (example: Chromosome 1) of distances to chromosome ends for TAG and WGD pairs, with centromere region and telomere threshold highlighted.

# MCSCANX

## Duplicates Prediction

To systematically classify duplicated genes, we used MCScanX, a widely adopted tool for detecting and categorizing gene duplications based on sequence similarity and chromosomal context. MCScanX first identifies homologous gene pairs through all-vs-all BLAST searches, then scans the genome for collinear blocks—regions where multiple homologous genes are arranged in the same order and orientation. Based on the arrangement and proximity of these homologs, MCScanX assigns each gene pair to one of several duplication categories:

- (Tandem) duplicates: Adjacent homologous genes on the same chromosome, with no or very few intervening genes.
- (Proximal) duplicates: Homologous genes located on the same chromosome but separated by a small number of non-homologous genes.
- (Dispersed) duplicates: Homologous genes located on different chromosomes or far apart on the same chromosome, not fitting other categories.
- (Segmental/WGD) duplicates: Homologous gene pairs that are part of larger collinear blocks, typically resulting from whole-genome or large-scale segmental duplications.

This classification enables a comprehensive view of the duplication landscape, distinguishing between local (tandem/proximal) and large-scale (segmental/WGD) events, and providing a robust framework for downstream evolutionary and functional analyses.

## Synteny

MCScanX also detects collinear (syntenic) blocks—conserved regions of gene order between different chromosomal segments—by chaining together homologous gene pairs that are co-linear and co-oriented. These blocks represent the remnants of ancient duplication events, such as whole-genome duplications, and are key to understanding genome evolution.

To visualize the extent and distribution of these collinear blocks, we used the MCScanX output to generate synteny plots. In these plots, each block is represented as a set of connected gene pairs spanning two chromosomal regions, highlighting conserved gene order and orientation. The resulting synteny maps provide a genome-wide view of duplicated segments, revealing patterns of large-scale structural conservation and rearrangement. These visualizations are instrumental in interpreting the evolutionary history of the soybean genome and in identifying regions of functional or evolutionary significance.

# Optimizations 

Processing over 2.4 million gene pairs required substantial workflow optimization to ensure both speed and reliability. We implemented parallelization throughout the pipelines, distributing BLAST searches, alignments, and Ks calculations across multiple CPU cores to dramatically reduce runtime. Batch processing was used to split large tasks into manageable chunks, allowing for efficient memory usage and easier error recovery. Checkpointing was integrated at key stages, so that intermediate results could be saved and the workflow could resume from the last successful step in case of interruption or failure. Additionally, all command-line steps were carefully optimized for performance, including the use of efficient file formats, streamlined data parsing, and robust error handling. These optimizations enabled us to process large genomic datasets reproducibly and efficiently, and to easily adapt the workflow for new species or parameter settings.

# References

Duan, X., Zhang, K., Duanmu, H., & Yu, Y. (2023). The myosin family genes in soybean: Genome-wide identification and expression analysis. South African Journal of Botany, 160, 338–346. https://doi.org/10.1016/j.sajb.2023.06.054
Miura, K., Toh, H., Hirakawa, H., Sugii, M., Murata, M., Nakai, K., Tashiro, K., Kuhara, S., Azuma, Y., & Shirai, M. (2008). Genome-wide analysis of Chlamydophila pneumoniae gene expression at the late stage of infection. DNA Research, 15(2), 93–102. https://doi.org/10.1093/dnares/dsn001
Kim, K. D., El Baidouri, M., Abernathy, B., Iwata-Otsubo, A., Chavarro, C., Gonzales, M., … Jackson, S. A. (2015). A comparative epigenomic analysis of polyploidy-derived genes in soybean and common bean. Plant Physiology, 168(4), 1433–1447. https://doi.org/10.1104/pp.15.00408
Yang, Y., Wang, J., & Di, J. (2013). Comparative inference of duplicated genes produced by polyploidization in the soybean genome. International Journal of Genomics, 2013, Article 810403. https://doi.org/10.1155/2013/810403
Roulin, A., Auer, P. L., Libault, M., Schlueter, J., Farmer, A., May, G., Stacey, G., Doerge, R. W., & Jackson, S. A. (2013). The fate of duplicated genes in a polyploid plant genome. The Plant Journal, 73(1), 143–153. https://doi.org/10.1111/tpj.12026
Schmutz, J., Cannon, S. B., Schlueter, J., Ma, J., Mitros, T., Nelson, W., Hyten, D. L., Song, Q., Thelen, J. J., Cheng, J., Xu, D., Hellsten, U., May, G. D., Yu, Y., Sakurai, T., Umezawa, T., Bhattacharyya, M. K., Sandhu, D., … Jackson, S. A. (2010). Genome sequence of the palaeopolyploid soybean. Nature, 463, 178–183. https://doi.org/10.1038/nature08670