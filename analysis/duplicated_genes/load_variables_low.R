

# -- data prep
families_df<-read.table('../../output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_evalue1e-10_wcol12_network.tsv', header=TRUE)
dup_df<-read.csv('../../output/info/duplicated_genes_info_id30_qcov50_scov50_evalue1e-10_wcol12.csv')
prot_df<-read.csv('../../output/info/protein_info_longest.csv')

dup_full_df <- dup_df[!is.na(dup_df$chromosome) & !is.na(dup_df$start_pos) & !is.na(dup_df$end_pos), ]

tags_df<-read.table('../../output/duplication_classes/TAGS/low/TAGs_1.tsv', header=TRUE)
singletons_df<-read.csv('../../output/info/singletons_genes_info_low.csv',header=T)


# ---------------- plots invariants ------------------

base_theme <- theme_minimal(base_size = 14)

custom_theme <- theme(
  plot.title = element_text( size = 18, hjust = 0.5),
  axis.title = element_text( size = 14),
  axis.text = element_text(size = 12),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "right"
)

# -- chormosome colors
chrom_colors <- c(
  "#1f77b4", "#1ca2c7", "#17becf", "#2ca02c", "#33cc33",
  "#98df8a", "#bcbd22", "#dbdb8d", "#ffbb78", "#ff7f0e",
  "#ff6600", "#d62728", "#ff9896", "#c5b0d5", "#9467bd",
  "#8c564b", "#c49c94", "#7f7f7f", "#aec7e8", "#393b79"
)

