# ------------------------------------------------------------------
# Script Name: 13_analyze_ks_vs_te_distance.R
#
# Description: 
#   Analyzes the relationship between Gene Age (Ks) and TE Proximity.
#   - GOAL: Test if older genes (High Ks) are located further from TEs 
#           than younger genes (Low Ks).
#   - INPUT: Gene Pairs (Ks), Gene Coordinates, TE GFF3.
#   - VIZ: Scatter Plot with Trend Lines.
#
# Input:  
#   1. ../data/TEAnnotationFinal.gff3
#   2. ../data/gene_pairs_classified.tsv
#   3. ../data/protein_info_longest.csv
# Output: ../output/TE_ks_vs_distance_scatter.png
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
print("Loading TEs...")
te_df <- read.table(te_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)

# NO FILTERING: All rows are included
te_all <- te_df %>% dplyr::select(V1, V4, V5) 
colnames(te_all) <- c("Chr", "Start", "End")

print(paste("Filtering stats: Keeping all", nrow(te_all), "TEs (classified and unclassified) for distance calculation."))

te_gr <- GRanges(seqnames = te_all$Chr, ranges = IRanges(start = te_all$Start, end = te_all$End))

# ==================================================================
# 2. LOAD GENE PAIRS (With Ks)
# ==================================================================
pairs_path  <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/gene_pairs_classified.tsv"
print("Loading Gene Pairs...")
pairs_df <- read.table(pairs_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# ==================================================================
# 3. LOAD GENE COORDINATES
# ==================================================================
coords_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/protein_info_longest.csv"
print("Loading Gene Coordinates...")
coords_df <- read.csv(coords_path, stringsAsFactors=FALSE) %>%
  dplyr::select(GeneID = peptide_id, Chr = chromosome, Start = start_pos, End = end_pos) %>%
  filter(!is.na(Start))

# ==================================================================
# 4. MERGE DATA (Flatten Pairs to Genes)
# ==================================================================
# A pair has 2 genes. We want to analyze them individually but keep the Pair's Ks.
gene1_data <- pairs_df %>% dplyr::select(GeneID = gene1, Ks = ks, Type = duplication_type)
gene2_data <- pairs_df %>% dplyr::select(GeneID = gene2, Ks = ks, Type = duplication_type)

genes_with_ks <- bind_rows(gene1_data, gene2_data) %>%
  inner_join(coords_df, by = "GeneID") %>%
  # Filter out scaffolds to match TE data
  filter(as.character(Chr) %in% unique(as.character(te_df$V1)))

# ==================================================================
# 5. CALCULATE DISTANCE TO NEAREST TE
# ==================================================================
print("Calculating distances...")
gene_gr <- GRanges(
  seqnames = genes_with_ks$Chr,
  ranges = IRanges(start = genes_with_ks$Start, end = genes_with_ks$End)
)

nearest_hits <- distanceToNearest(gene_gr, te_gr, ignore.strand=TRUE)
genes_with_ks$Dist_to_TE <- mcols(nearest_hits)$distance

# ==================================================================
# 6. VISUALIZATION
# ==================================================================
# Filter unrealistic Ks values (e.g., > 5.0 is usually saturation noise)
plot_data <- genes_with_ks %>% 
  filter(Ks < 3.0)

# Colors for types 
type_colors <- c("WGD" = "#7A0177",        
                 "TAG" = "#084594",        
                 "Dispersed" = "#FA9FB5")  

p <- ggplot(plot_data, aes(x = Ks, y = Dist_to_TE, color = Type)) +
  
  # A. Scatter Points (Small and transparent to see density)
  geom_point(alpha = 0.3, size = 0.8) +
  
  # B. Smooth Trend Line (GAM) to show the pattern
  geom_smooth(method = "gam", color = "black", size = 0.8) +
  
  # Facet by Type to see distinct evolutionary patterns
  facet_wrap(~Type, ncol = 1) +
  
  scale_color_manual(values = type_colors) +
  
  # Log scale for Distance (Y-axis)
  scale_y_continuous(trans = "log1p", 
                     breaks = c(0, 100, 1000, 10000, 50000),
                     labels = c("0", "100bp", "1kb", "10kb", "50kb")) +
  
  theme_bw() +
  
  labs(title = "Evolutionary Dynamics: Gene Age vs. TE Proximity",
       subtitle = "Distance to nearest ALL TEs (Classified and Unclassified)",
       x = "Synonymous Substitution Rate (Ks) ~ Age",
       y = "Distance to Nearest TE (bp)") +
  
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

# Save
output_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_ks_vs_distance_scatter.png"
ggsave(output_path, plot = p, width = 8, height = 10, dpi = 300)

print(paste("Plot saved to:", output_path))