# ------------------------------------------------------------------
# Script Name: 11_analyze_nesting_patterns.R
#
# Description: 
#   Determines specific Nesting Relationships ("Who is inside Whom?").
#   - METHOD: Pairwise Overlap Calculation using GenomicRanges.
#   - METRIC: % of Family A (Row) that physically overlaps Family B (Col).
#   - SCOPE: Top 20 most abundant families (to keep matrix readable).
#   - OUTPUT: Heatmap of Nesting Percentages.
#
# Input:  /data/TEAnnotationFinal.gff3
# Output: /output/TE_nesting_heatmap.png
# ------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ==================================================================
# 1. LOAD & CLEAN DATA
# ==================================================================
input_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/TEAnnotationFinal.gff3"

gff_data <- read.table(input_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)
raw_data <- gff_data %>% dplyr::select(V1, V3, V4, V5)
colnames(raw_data) <- c("Chr", "Full_Name", "Start", "End")

# Clean and Parse
clean_data <- raw_data %>%
  tidyr::separate(Full_Name, into = c("Class_Raw", "Order", "Superfamily"), sep = "/", fill = "right") %>%
  mutate(
    Superfamily = ifelse(Superfamily == "-" | is.na(Superfamily) | Superfamily == "", "Unknown", Superfamily),
    # Create Display Name
    DisplayName = paste0(Order, "/", Superfamily)
  ) %>%
  filter(!grepl("scaffold|contig|chl|mito|Unknown", Chr, ignore.case = TRUE)) %>%
  filter(Superfamily != "Unknown")

# ==================================================================
# 2. DEFINE TOP FAMILIES (For Computation Speed & readability)
# ==================================================================
top_families <- clean_data %>%
  count(DisplayName) %>%
  arrange(desc(n)) %>%
  slice_head(n = 20) %>%
  pull(DisplayName)

print("Analyzing nesting for Top 20 Superfamilies:")
print(top_families)

target_data <- clean_data %>% filter(DisplayName %in% top_families)

# ==================================================================
# 3. CALCULATE PAIRWISE OVERLAPS
# ==================================================================
# Convert to GRanges list (one GRanges object per Family)
gr_list <- split(
  GRanges(seqnames = target_data$Chr, ranges = IRanges(target_data$Start, target_data$End)),
  target_data$DisplayName
)

results_list <- list()
counter <- 1

# Loop through every pair
for (guest_name in names(gr_list)) {
  
  guest_gr <- gr_list[[guest_name]]
  total_guest_bp <- sum(width(guest_gr))
  
  for (host_name in names(gr_list)) {
    
    host_gr <- gr_list[[host_name]]
    
    # Calculate Intersection
    # findOverlaps gets the indices of overlaps
    hits <- findOverlaps(guest_gr, host_gr)
    
    # If looking at Self-Nesting (Gypsy vs Gypsy), ignore "Self-Hits" (Identity)
    if (guest_name == host_name) {
      hits <- hits[queryHits(hits) != subjectHits(hits)]
    }
    
    # Calculate the physical BP of the overlap
    if (length(hits) > 0) {
      overlaps <- pintersect(guest_gr[queryHits(hits)], host_gr[subjectHits(hits)])
      overlap_bp <- sum(width(overlaps))
    } else {
      overlap_bp <- 0
    }
    
    # Calculate % of Guest that is inside Host
    pct_overlap <- (overlap_bp / total_guest_bp) * 100
    
    results_list[[counter]] <- data.frame(
      Guest = guest_name,
      Host = host_name,
      Overlap_Pct = pct_overlap
    )
    counter <- counter + 1
  }
}

matrix_df <- do.call(rbind, results_list)

# ==================================================================
# 4. VISUALIZATION
# ==================================================================
# Sort axes by abundance (using the top_families order)
matrix_df$Guest <- factor(matrix_df$Guest, levels = rev(top_families)) # Rev for Y-axis top-down
matrix_df$Host <- factor(matrix_df$Host, levels = top_families)

p <- ggplot(matrix_df, aes(x = Host, y = Guest, fill = Overlap_Pct)) +
  
  geom_tile(color = "white") +
  
  # Color: White (0%) -> Light Blue -> Dark Blue -> Purple (High %)
  scale_fill_gradientn(
    colors = c("white", "#9ECAE1", "#084594", "#7A0177"),
    values = rescale(c(0, 1, 5, 20, 50)), # Bias color scale to show low % clearly
    name = "% of Guest\ninside Host"
  ) +
  
  theme_bw() +
  
  labs(title = "TE Nesting Matrix (Who is inside Whom?)",
       subtitle = "Color shows the percentage of the 'Guest' (Row) that physically overlaps the 'Host' (Column).",
       x = "HOST (Container)",
       y = "GUEST (Nested Element)") +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face="bold"),
    axis.text.y = element_text(size = 9, face="bold"),
    panel.grid = element_blank()
  )

# ==================================================================
# 5. SAVE OUTPUT
# ==================================================================
output_file <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_nesting_heatmap.png"
ggsave(output_file, plot = p, width = 12, height = 10, dpi = 300)

print(paste("Nesting Heatmap saved to:", output_file))