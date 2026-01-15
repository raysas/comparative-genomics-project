# ------------------------------------------------------------------
# Script Name: 12_analyze_gene_te_proximity.R
#
# Description: 
#   Integrates TE data with Gene Duplication Classification.
#   - GOAL: Determine if TEs are differentially distributed around 
#           WGD, Tandem, Dispersed, and Singleton genes.
#   - LOGIC:
#       1. Parse 'gene_pairs_classified.tsv' to get types for paired genes.
#       2. Parse 'protein_info_longest.csv' to get coordinates for ALL genes.
#       3. Identify "Singletons" (Genes in coord file but NOT in pairs file).
#       4. Calculate distance to nearest ALL TEs.
#   - OUTPUT: Boxplot and Stats CSV.
#
# Input:  
#   1. ../data/TEAnnotationFinal.gff3
#   2. ../data/gene_pairs_classified.tsv  (gene1, gene2, duplication_type)
#   3. ../data/protein_info_longest.csv   (peptide_id, chromosome, start_pos, end_pos)
# Output: 
#   - ../output/TE_gene_proximity_boxplot.png
#   - ../output/TE_gene_proximity_stats.csv
# ------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")

library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(scales)

# ==================================================================
# 1. LOAD TE DATA (INCLUDING UNKNOWNS)
# ==================================================================
te_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/TEAnnotationFinal.gff3"

print("Loading ALL TEs (including Unknowns)...")
# Read raw table 
te_raw <- read.table(te_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)

# NO FILTERING: All rows are included
te_all <- te_raw %>% dplyr::select(V1, V4, V5) 
colnames(te_all) <- c("Chr", "Start", "End")

# Convert to GRanges
te_gr <- GRanges(
  seqnames = te_all$Chr, 
  ranges = IRanges(start = te_all$Start, end = te_all$End)
)

# ==================================================================
# 2. LOAD & PROCESS GENE DATA
# ==================================================================
pairs_path  <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/gene_pairs_classified.tsv"
coords_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/protein_info_longest.csv"

# A. Load Coordinates (The Master List of Genes)
print("Loading Gene Coordinates...")
coords_df <- read.csv(coords_path, stringsAsFactors=FALSE) %>%
  dplyr::select(peptide_id, chromosome, start_pos, end_pos) %>%
  dplyr::rename(GeneID = peptide_id, Chr = chromosome, Start = start_pos, End = end_pos) %>%
  filter(!is.na(Start) & !is.na(End)) # Safety check

# B. Load Classifications (The Pairs)
print("Loading Gene Pairs...")
pairs_df <- read.table(pairs_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# C. Reshape Pairs to Gene List (Long Format)
classified_genes <- bind_rows(
  pairs_df %>% dplyr::select(GeneID = gene1, Type = duplication_type),
  pairs_df %>% dplyr::select(GeneID = gene2, Type = duplication_type)
) %>%
  distinct(GeneID, .keep_all = TRUE)

# D. Merge and Identify Singletons
print("Merging and classifying Singletons...")
gene_data <- left_join(coords_df, classified_genes, by="GeneID") %>%
  mutate(
    Type = ifelse(is.na(Type), "Singleton", Type) # Identify Singletons
  ) %>%
  # Filter out scaffolds not present in TE file
  filter(as.character(Chr) %in% unique(as.character(te_all$Chr)))

print(paste("Total Genes Analyzed:", nrow(gene_data)))
print(table(gene_data$Type))

# Create Gene GRanges
gene_gr <- GRanges(
  seqnames = gene_data$Chr,
  ranges = IRanges(start = gene_data$Start, end = gene_data$End),
  gene_id = gene_data$GeneID,
  type = gene_data$Type
)

# ==================================================================
# 3. CALCULATE SPATIAL PROXIMITY
# ==================================================================
print("Calculating distance to nearest ALL TEs...")

# METHOD: distanceToNearest()
nearest_hits <- distanceToNearest(gene_gr, te_gr, ignore.strand=TRUE)

# Extract distances and add to dataframe
gene_data$Dist_to_TE <- mcols(nearest_hits)$distance

# ==================================================================
# 4. STATISTICAL SUMMARY
# ==================================================================
summary_stats <- gene_data %>%
  group_by(Type) %>%
  summarise(
    Count = n(),
    Median_Dist_bp = median(Dist_to_TE),
    Mean_Dist_bp = mean(Dist_to_TE),
    Q3_Dist_bp = quantile(Dist_to_TE, 0.75),
    Overlap_Percent = round(sum(Dist_to_TE == 0) / n() * 100, 2)
  ) %>%
  arrange(Median_Dist_bp)

print("--- Summary Statistics ---")
print(summary_stats)

# Save Stats
write.csv(summary_stats, "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_gene_proximity_stats_all_TEs.csv", row.names=FALSE)

# Statistical Test (Pairwise Wilcoxon)
stat_test <- pairwise.wilcox.test(gene_data$Dist_to_TE, gene_data$Type, p.adjust.method = "bonferroni")
print("Significance Test (P-values):")
print(stat_test)

# ==================================================================
# 5. VISUALIZATION
# ==================================================================
# Reorder Factor by Median Distance for clean plotting
gene_data$Type <- factor(gene_data$Type, levels = summary_stats$Type)

# Define Colors for Gene Types (Custom Colors)
type_colors <- c("WGD" = "#7A0177",       
                 "Singleton" = "#084594",  
                 "TAG" = "#FDD0A2",        
                 "Dispersed" = "#A6BDDB")  

p <- ggplot(gene_data, aes(x = Type, y = Dist_to_TE, fill = Type)) +
  
  # A. The Boxplot
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  
  # B. Log Scale Y-Axis
  scale_y_continuous(trans = "log1p", 
                     breaks = c(0, 10, 100, 1000, 10000, 50000, 100000),
                     labels = c("0 (Overlap)", "10 bp", "100 bp", "1 kb", "10 kb", "50 kb", "100 kb")) +
  
  scale_fill_manual(values = type_colors) +
  
  theme_bw() +
  
  labs(title = "Gene Proximity to All TEs",
       subtitle = "Distance (Edge-to-Edge) to nearest TE (including Unclassified/Unknown TEs).",
       x = "Duplication Type",
       y = "Distance to Nearest TE (bp)") +
  
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10),
        plot.title = element_text(face="bold"),
        panel.grid.major.x = element_blank())

# Save Plot
out_png <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_gene_proximity_boxplot_all_TEs.png"
ggsave(out_png, plot = p, width = 8, height = 7, dpi = 300)

print(paste("Plot saved to:", out_png))