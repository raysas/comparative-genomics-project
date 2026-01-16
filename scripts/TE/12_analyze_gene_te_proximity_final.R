# ------------------------------------------------------------------
# Script Name: 12_analyze_gene_te_proximity_final.R
# Description: Merges duplication lists and calculates TE proximity
#              for both HIGH and LOW stringencies with stats on plots.
# ------------------------------------------------------------------

# 1. INSTALL/LOAD NECESSARY PACKAGES
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
if (!require("ggpubr", quietly = TRUE)) install.packages("ggpubr")

library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(scales)
library(ggpubr)

# Paths - Adjust these if your folder structure is different
data_dir <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/"
out_dir  <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/"
te_path  <- paste0(data_dir, "TEAnnotationFinal.gff3")
coords_path <- paste0(data_dir, "protein_info_longest.csv")

# ==================================================================
# 1. LOAD TE DATA (COMMON)
# ==================================================================
print("Loading TE data...")
te_raw <- read.table(te_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)
te_all <- te_raw %>% dplyr::select(V1, V4, V5) 
colnames(te_all) <- c("Chr", "Start", "End")

te_gr <- GRanges(
  seqnames = te_all$Chr, 
  ranges = IRanges(start = te_all$Start, end = te_all$End)
)

# Load Coordinates
coords <- read.csv(coords_path) %>%
  dplyr::select(peptide_id, chromosome, start_pos, end_pos) %>%
  dplyr::rename(GeneID = peptide_id, Chr = chromosome, Start = start_pos, End = end_pos)

# ==================================================================
# 2. RUN ANALYSIS FOR BOTH STRINGENCIES
# ==================================================================
stringencies <- c("high", "low")

for (strig in stringencies) {
  print(paste("--- PROCESSING:", toupper(strig), "STRINGENCY ---"))
  
  # A. Merge individual lists into one classification table
  sings <- read.table(paste0(data_dir, "singletons_", strig, ".txt"), stringsAsFactors=F)$V1
  tags  <- read.table(paste0(data_dir, "TAGs_", strig, ".txt"), stringsAsFactors=F)$V1
  dups  <- read.table(paste0(data_dir, "duplicated_genes_", strig, ".txt"), stringsAsFactors=F)$V1
  
  class_table <- bind_rows(
    data.frame(GeneID = sings, Type = "Singleton"),
    data.frame(GeneID = tags,  Type = "TAG"),
    data.frame(GeneID = dups,  Type = "Duplicated")
  )
  
  gene_data <- inner_join(coords, class_table, by = "GeneID") %>%
    filter(as.character(Chr) %in% unique(as.character(te_all$Chr)))

  # B. Calculate Distance
  gene_gr <- GRanges(
    seqnames = gene_data$Chr,
    ranges = IRanges(start = gene_data$Start, end = gene_data$End),
    type = gene_data$Type
  )
  
  print("Calculating distances...")
  nearest_hits <- distanceToNearest(gene_gr, te_gr, ignore.strand=TRUE)
  gene_data$Dist_to_TE <- mcols(nearest_hits)$distance

  # C. Summary Stats
  summary_stats <- gene_data %>%
    group_by(Type) %>%
    summarise(Median = median(Dist_to_TE), Count = n()) %>%
    arrange(Median)
  
  write.csv(summary_stats, paste0(out_dir, "proximity_stats_", strig, ".csv"), row.names=F)

  # D. Plotting with Stats
  gene_data$Type <- factor(gene_data$Type, levels = summary_stats$Type)
  
  type_colors <- c("Duplicated" = "#7A0177", "Singleton" = "#084594", "TAG" = "#FDD0A2")
  my_comparisons <- list( c("Duplicated", "Singleton"), c("Duplicated", "TAG"), c("Singleton", "TAG") )

  p <- ggplot(gene_data, aes(x = Type, y = Dist_to_TE, fill = Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
    scale_y_continuous(trans = "log1p", 
                       breaks = c(0, 100, 1000, 10000, 100000),
                       labels = c("0", "100 bp", "1 kb", "10 kb", "100 kb")) +
    scale_fill_manual(values = type_colors) +
    theme_bw() +
    labs(title = paste("Gene-TE Proximity (", toupper(strig), ")"),
         subtitle = "Wilcoxon test: ***p<0.001, **p<0.01, *p<0.05, ns=non-significant",
         x = "Duplication Type", y = "Distance to Nearest TE (bp)") +
    theme(legend.position = "none", axis.text.x = element_text(size=12, face="bold"))

  ggsave(paste0(out_dir, "TE_proximity_plot_", strig, ".png"), plot = p, width = 8, height = 7)
  print(paste("Completed", strig))
}