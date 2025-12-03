BiocManager::install('clusterProfiler')
BiocManager::install('org.Gmax.eg.db')
library(clusterProfiler)
library(org.Gmax.eg.db)

largest_family=readLines('../../output/gene_lists/largest_family_id50_cov70_evalue1e-10',n=1)

# -- running clusterProfiler enrichment analysis
gene.df <- bitr(largest_family, fromType = "Gmax_Gene_ID", 
                toType = c("ENTREZID"),
                OrgDb = org.Gmax.eg.db)
ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Gmax.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)