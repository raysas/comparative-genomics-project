# longest_family_members
# # -- save it in text file
# 
# 
# # -- distribution of family members across chromosomes
# longest_family_dup_df<-dup_full_df[dup_full_df$peptide_id %in% longest_family_members, ]
# table(longest_family_dup_df$chromosome)
# 
# # -- how much of them are TAGs 
# longest_family_tags_df<-tags_df[tags_df$peptide_id %in% longest_family_members, ]
# table(longest_family_tags_df$TAG)
# length(table(longest_family_tags_df$TAG))
# order(table(longest_family_tags_df$TAG), decreasing=TRUE)

# -- most of them are not (713 out of total)
# -- the rest are split to 168 different TAG blocks

edgelist<-read.table('../../output/similarity_edgelists/filtered_blast_results_id50_qcov70_scov70_evalue1e-10_wcol12_network.tsv')
largest_fam<-read.table('../../output/gene_lists/largest_family/largest_family_high.txt', header=FALSE)
largest_family_edgelist<-edgelist[edgelist$V1 %in% largest_fam$V1 & edgelist$V2 %in% largest_fam$V1, ]
library(igraph)
g<-graph_from_data_frame(largest_family_edgelist, directed=FALSE)
# -- large graph can not be visualized

# -- viz
plot(g, vertex.size=5, vertex.label=NA)
  
# -- save graph in file
library(ggraph)
library(ggplot2)
p <- ggraph(g, layout = 'fr') +
  geom_edge_link(alpha = 0.5) +
  geom_node_point(size = 3, color = 'blue') +
  theme_void()

write_graph(g, file = "data/largest_cluster.graphml", format = "graphml")
