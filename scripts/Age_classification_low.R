library(dplyr)
library(readr)


df <- read_delim(
  "../output/ks_results_low/ks_with_GLYMA_ID_GO_low.tsv",
  delim = "\t",
  col_types = cols()
)

# EXTRAIRE les gènes Old 
old_genes <- df %>%
  filter(age_class == "Old") %>%
  pull(gene) %>%
  unique()

# EXTRAIRE les gènes Recent
recent_genes <- df %>%
  filter(age_class == "Recent") %>%
  pull(gene) %>%
  unique()


write.table(old_genes,
            "../output/ks_results_low/old_gene_list_low.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(recent_genes,
            "../output/ks_results_low/recent_gene_list_low.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
