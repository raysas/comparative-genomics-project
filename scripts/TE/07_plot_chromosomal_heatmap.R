# ------------------------------------------------------------------
# Script Name: 07_plot_chromosomal_heatmap.R
#
# Description: 
#   Visualizes the Spatial Distribution of TEs across the genome (Heatmap).
#   - USES EXACT CHROMOSOME LENGTHS from reference file.
#   - METHOD: Uses `GenomicRanges` to calculate TRUE COVERAGE.
#   - METRIC: Coverage % in 100kb windows.
#   - VIZ: Heatmap (White -> Blue -> Purple).
#   - INCLUDES: Class I, Class II, and Unknowns.
#
# Input:  /data/TEAnnotationFinal.gff3
#         /data/chr_lengths.tsv
# Output: /output/TE_chromosomal_heatmap.png
# ------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(GenomicRanges) 

# ==================================================================
# 1. DATA LOADING & PRE-PROCESSING
# ==================================================================
gff_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/TEAnnotationFinal.gff3"
len_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/chr_lengths.tsv"

# A. Load True Chromosome Lengths
# Format: chromosome <tab> length
chr_ref <- read.table(len_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(chr_ref) <- c("Chr", "Total_Length")

# Convert to a named vector for GenomicRanges
seq_lengths <- setNames(chr_ref$Total_Length, as.character(chr_ref$Chr))

# B. Load TE Data
gff_data <- read.table(gff_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)
raw_data <- gff_data %>% select(V1, V3, V4, V5)
colnames(raw_data) <- c("Chr", "Full_Name", "Start", "End")

# Clean and Parse
clean_data <- raw_data %>%
  separate(Full_Name, into = c("Class_Raw", "Order", "Superfamily"), sep = "/", fill = "right") %>%
  mutate(
    Class_Panel = case_when(
      grepl("Class II", Class_Raw) ~ "Class II (DNA Transposons)",
      grepl("Class I", Class_Raw) ~ "Class I (Retrotransposons)",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!grepl("scaffold|contig|chl|mito|Unknown", Chr, ignore.case = TRUE)) %>%
  filter(Start > 0 & End > Start)


# ==================================================================
# 2. DEFINE GENOME WINDOWS (Using True Lengths)
# ==================================================================
BIN_SIZE <- 100000 # 100 kb 

# Create Tile Genome using the REFERENCE lengths
# This creates bins all the way to the true end of the chromosome
windows_gr <- tileGenome(seqlengths = seq_lengths, 
                         tilewidth = BIN_SIZE, 
                         cut.last.tile.in.chrom = TRUE)

# ==================================================================
# 3. CALCULATE TRUE COVERAGE
# ==================================================================
calc_class_coverage <- function(data_subset, windows, seq_lens) {
  
  # Create GRanges for TEs with explicit seqlengths
  te_gr <- GRanges(
    seqnames = factor(data_subset$Chr, levels = names(seq_lens)),
    ranges = IRanges(start = data_subset$Start, end = data_subset$End),
    seqlengths = seq_lens
  )
  
  # Calculate Coverage Vector (Rle)
  cov_rle <- coverage(te_gr)
  
  # Create Boolean Coverage (Union)
  cov_bool <- cov_rle > 0 
  
  # Calculate mean (proportion covered) per window
  res_gr <- binnedAverage(windows, cov_bool, "Density_Prop")
  
  # Convert to DF
  df <- as.data.frame(res_gr, row.names = NULL)
  df$Density_Pct <- df$Density_Prop * 100
  return(df)
}

# Run for Class I
c1_data <- clean_data %>% filter(Class_Panel == "Class I (Retrotransposons)")
c1_res <- calc_class_coverage(c1_data, windows_gr, seq_lengths)
c1_res$Class_Panel <- "Class I (Retrotransposons)"

# Run for Class II
c2_data <- clean_data %>% filter(Class_Panel == "Class II (DNA Transposons)")
c2_res <- calc_class_coverage(c2_data, windows_gr, seq_lengths)
c2_res$Class_Panel <- "Class II (DNA Transposons)"

# Run for Unknown
unk_data <- clean_data %>% filter(Class_Panel == "Unknown")
unk_res <- calc_class_coverage(unk_data, windows_gr, seq_lengths)
unk_res$Class_Panel <- "Unknown"

# Combine
final_df <- rbind(c1_res, c2_res, unk_res)

# ==================================================================
# 4. PREPARE FOR PLOTTING
# ==================================================================
final_df$Chr_Num <- as.numeric(as.character(final_df$seqnames))
final_df$Chr_Factor <- factor(final_df$Chr_Num, levels = rev(sort(unique(final_df$Chr_Num))))
final_df$Position_Mb <- final_df$start / 1000000

# ==================================================================
# 5. VISUALIZATION
# ==================================================================
p <- ggplot(final_df, aes(x = Position_Mb, y = Chr_Factor, fill = Density_Pct)) +
  
  geom_tile() +
  
  facet_wrap(~Class_Panel, ncol = 3) +
  
  # White -> Blue -> Purple
  scale_fill_gradientn(colors = c("white", "#9ECAE1", "#084594", "#7A0177"),
                       name = "TE Coverage (%)",
                       limits = c(0, 100)) +
  
  scale_x_continuous(expand = c(0, 0)) + 
  
  theme_bw() +
  
  labs(title = "Chromosomal Distribution of TEs",
       subtitle = "TE Coverage (100kb windows)",
       x = "Genomic Position (Mb)",
       y = "Chromosome") +
  
  theme(strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

# ==================================================================
# 6. SAVE OUTPUT
# ==================================================================
output_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_chromosomal_heatmap.png"
ggsave(output_path, plot = p, width = 16, height = 8, dpi = 300)

print(paste("Plot saved to:", output_path))