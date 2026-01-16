# ------------------------------------------------------------------
# Script Name: 12_map_gene_classifications.R
# Description: Merges genomic coordinates with gene duplication categories (Singletons, TAGs, Duplicates).
# Ensures categories are mutually exclusive by removing TAGs from the general duplicate pool.
# ------------------------------------------------------------------

library(dplyr)

# 1. Load genomic coordinates (longest isoforms only)
coords_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/protein_info_longest.csv"
coords <- read.csv(coords_path) %>%
  select(peptide_id, chromosome, start_pos, end_pos) %>%
  rename(GeneID = peptide_id, Chr = chromosome, Start = start_pos, End = end_pos)

# 2. Function to merge coordinates and resolve category overlaps
merge_stringency <- function(level, coords_df) {
  data_dir <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/"
  
  # Load raw gene lists
  sings     <- read.table(paste0(data_dir, "singletons_", level, ".txt"), stringsAsFactors=F)$V1
  tags      <- read.table(paste0(data_dir, "TAGs_", level, ".txt"), stringsAsFactors=F)$V1
  all_dups  <- read.table(paste0(data_dir, "duplicated_genes_", level, ".txt"), stringsAsFactors=F)$V1
  
  # Ensure mutual exclusivity: Remove TAGs from the duplicate list
  non_tags <- setdiff(all_dups, tags)
  
  # Create classification table
  class_table <- bind_rows(
    data.frame(GeneID = sings,    Type = "Singleton"),
    data.frame(GeneID = tags,     Type = "TAG"),
    data.frame(GeneID = non_tags, Type = "Non-TAG Duplicate")
  )
  
  # Merge with coordinates and remove any duplicates
  merged <- inner_join(coords_df, class_table, by = "GeneID") %>%
    distinct(GeneID, .keep_all = TRUE)
  
  return(merged)
}

# 3. Process and export High/Low stringency files
gene_data_high <- merge_stringency("high", coords)
gene_data_low  <- merge_stringency("low", coords)

write.csv(gene_data_high, "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/merged_genes_high.csv", row.names=FALSE)
write.csv(gene_data_low,  "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/merged_genes_low.csv", row.names=FALSE)

# 4. Print summary statistics
print("Files successfully created:")
print(table(gene_data_high$Type))