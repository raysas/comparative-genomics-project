
dup_df<-read.csv('../../output/statistics/duplicated_genes_info.csv')
prot_df<-read.csv('../../output/statistics/protein_info_longest.csv')
edgelist<-read.table('../../output/similarity_edgelists/filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv')
families_df<-read.table('../../output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv', header=TRUE)
tags_info_df<-read.table('../../output/statistics/TAGs_ratios.tsv', header=TRUE)
tags_df<-read.table('../../output/duplication_classes/TAGS/TAGs_1.tsv', header=TRUE)


colnames(dup_df)
dup_df$chromosome <- factor(dup_df$chromosome)
prot_df$chromosome <- factor(prot_df$chromosome)

dup_full_df <- dup_df[!is.na(dup_df$chromosome) & !is.na(dup_df$start_pos) & !is.na(dup_df$end_pos), ]

# -- sort by chromosome number
dup_df<-dup_df[order(as.numeric(as.character(dup_df$chromosome))), ]
prot_df<-prot_df[order(as.numeric(as.character(prot_df$chromosome))), ]

# -- chormosome colors
chrom_colors <- c(
  "#1f77b4", "#1ca2c7", "#17becf", "#2ca02c", "#33cc33",
  "#98df8a", "#bcbd22", "#dbdb8d", "#ffbb78", "#ff7f0e",
  "#ff6600", "#d62728", "#ff9896", "#c5b0d5", "#9467bd",
  "#8c564b", "#c49c94", "#7f7f7f", "#aec7e8", "#393b79"
)
         
# ---------------- for largesdt family analysis -----------------
# -- get the longest family
count_df<-as.data.frame(table(families_df$family)) 
colnames(count_df)<-c('family','count')
count_df <-count_df[order(-count_df$count), ]
longest_family<-as.character(count_df$family[1])
longest_family_members<-families_df[families_df$family==longest_family, ]$geneName

# -- the the cluster of the longest family
longest_family_edgelist<-edgelist[edgelist$V1 %in% longest_family_members & edgelist$V2 %in% longest_family_members, ]
library(igraph)
g<-graph_from_data_frame(longest_family_edgelist, directed=FALSE)
# -- large graph can not be visualized