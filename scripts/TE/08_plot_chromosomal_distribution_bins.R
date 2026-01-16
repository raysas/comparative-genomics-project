# ------------------------------------------------------------------
# Script Name: 08_plot_chromosomal_distribution_bins.R
#
# Description: 
#   Visualizes TE distribution per chromosome at MULTIPLE RESOLUTIONS.
#   - LOOPS through bin sizes: 100kb, 200kb, 500kb, 1Mb, 2Mb, 5Mb.
#   - METRIC: Genomic Coverage (Mb) per bin.
#   - SCALING: Recalculates Y-axis limits dynamically for each resolution.
#   - OUTPUTS: A PDF (detailed) and PNG (grid) for EACH bin size.
#
# Input:  /data/TEAnnotationFinal.gff3
#         /data/chr_lengths.tsv
# Output: /output/TE_bars_{size}.pdf
#         /output/TE_bars_{size}.png
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

# ==================================================================
# 1. SETUP & DATA LOADING
# ==================================================================
gff_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/TEAnnotationFinal.gff3"
len_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/chr_lengths.tsv"

# Load Chromosome Lengths
chr_lengths_df <- read.table(len_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(chr_lengths_df) <- c("Chr_Num", "Total_Length_bp") 

# Load GFF3
gff_data <- read.table(gff_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)
raw_data <- gff_data %>% select(V1, V3, V4, V5, V9)
colnames(raw_data) <- c("Chr", "Full_Name", "Start", "End", "Attributes")

# Clean and Parse Data (Done once)
clean_data <- raw_data %>%
  separate(Full_Name, into = c("Class_Raw", "Order", "Superfamily"), sep = "/", fill = "right") %>%
  mutate(
    Class_Group = case_when(
      grepl("Class II", Class_Raw) ~ "Class II",
      grepl("Class I", Class_Raw) ~ "Class I",
      TRUE ~ "Unknown"
    ),
    Length_Str = str_extract(Attributes, "Length=[0-9]+"), 
    Length = as.numeric(gsub("Length=", "", Length_Str)),
    Length = ifelse(is.na(Length), End - Start + 1, Length),
    Order = ifelse(Order == "-" | is.na(Order) | Order == "", "Unknown", Order)
  ) %>%
  filter(!grepl("scaffold|contig|chl|mito|Unknown", Chr, ignore.case = TRUE)) 

clean_data$Chr_Num <- as.numeric(clean_data$Chr)
clean_data <- clean_data %>% filter(!is.na(Chr_Num)) 

# ==================================================================
# 2. DEFINE BIN SIZES TO COMPARE
# ==================================================================
# The script will run the analysis for every size in this list
bin_configs <- list(
  "100kb" = 100000,
  "200kb" = 200000,
  "500kb" = 500000,
  "1Mb"   = 1000000,
  "2Mb"   = 2000000,
  "5Mb"   = 5000000
)

# ==================================================================
# 3. COLOR GENERATION (Consistent across all plots)
# ==================================================================
# Calculate global abundance once to lock colors
global_order_abundance <- clean_data %>%
  group_by(Class_Group, Order) %>%
  summarise(Total = sum(Length), .groups='drop') %>%
  arrange(Class_Group, desc(Total))

c1_orders <- global_order_abundance %>% filter(Class_Group == "Class I") %>% pull(Order)
c2_orders <- global_order_abundance %>% filter(Class_Group == "Class II") %>% pull(Order)
unk_orders <- global_order_abundance %>% filter(Class_Group == "Unknown") %>% pull(Order)

all_orders <- c(c1_orders, c2_orders, unk_orders)
clean_data$Order <- factor(clean_data$Order, levels = all_orders)

get_gradient <- function(colors_vec, n) {
  if (n == 0) return(c())
  colorRampPalette(colors_vec)(n)
}

cols_c1 <- get_gradient(c("#7A0177", "#FA9FB5"), length(c1_orders))
if(length(cols_c1)>0) names(cols_c1) <- c1_orders
cols_c2 <- get_gradient(c("#084594", "#9ECAE1"), length(c2_orders))
if(length(cols_c2)>0) names(cols_c2) <- c2_orders
cols_unk <- rep("#999999", length(unk_orders))
if(length(cols_unk)>0) names(cols_unk) <- unk_orders

te_colors <- c(cols_c1, cols_c2, cols_unk)

# ==================================================================
# 4. MASTER LOOP (Iterate through Bin Sizes)
# ==================================================================

for (bin_name in names(bin_configs)) {
  
  CURRENT_BIN_SIZE <- bin_configs[[bin_name]]
  print(paste("Processing Bin Size:", bin_name, "...", CURRENT_BIN_SIZE, "bp"))
  
  # --- A. Binning Logic ---
  binned_data <- clean_data %>%
    mutate(
      Bin_Start = floor(Start / CURRENT_BIN_SIZE) * CURRENT_BIN_SIZE,
      Bin_Label = Bin_Start / 1000000 # Keep axis in Mb
    ) %>%
    group_by(Chr_Num, Bin_Label, Class_Group, Order) %>%
    summarise(Bin_Coverage_bp = sum(Length), .groups = 'drop') %>%
    mutate(Bin_Coverage_Mb = Bin_Coverage_bp / 1000000)
  
  # Add Labels for Faceting
  binned_data$Chr_Label <- paste("Chr", binned_data$Chr_Num)
  binned_data$Chr_Label <- factor(binned_data$Chr_Label, 
                                  levels = paste("Chr", sort(unique(binned_data$Chr_Num))))
  
  # --- B. Calculate Y-Axis Limit Specific to this Bin Size ---
  # Larger bins hold more TEs, so the Y-limit must scale up.
  total_bin_heights <- binned_data %>%
    group_by(Chr_Num, Bin_Label) %>%
    summarise(Total_Height = sum(Bin_Coverage_Mb), .groups = 'drop')
  
  GLOBAL_MAX_Y <- max(total_bin_heights$Total_Height) * 1.05
  
  # --- C. Output 1: Multi-Page PDF ---
  pdf_file <- paste0("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_bars_", bin_name, ".pdf")
  pdf(pdf_file, width = 12, height = 6)
  
  chromosomes <- sort(unique(binned_data$Chr_Num))
  for (chr in chromosomes) {
    chr_data <- binned_data %>% filter(Chr_Num == chr)
    
    # Get exact chr length
    chr_limit_mb <- chr_lengths_df$Total_Length_bp[chr_lengths_df$Chr_Num == chr] / 1000000
    
    p <- ggplot(chr_data, aes(x = Bin_Label, y = Bin_Coverage_Mb, fill = Order)) +
      geom_bar(stat = "identity", width = (CURRENT_BIN_SIZE/1000000) * 0.9) + 
      scale_fill_manual(values = te_colors) +
      scale_y_continuous(limits = c(0, GLOBAL_MAX_Y), expand = c(0, 0)) +
      scale_x_continuous(limits = c(0, chr_limit_mb), breaks = pretty_breaks(n=10), expand = c(0, 0)) +
      
      labs(title = paste0("TE Distribution on Chromosome ", chr, " (", bin_name, " bins)"),
           subtitle = paste0("Total Coverage (Mb) per ", bin_name, " Window"),
           x = "Position (Mb)", y = "Coverage (Mb)") +
      
      theme_bw() +
      theme(panel.grid = element_blank(), strip.text = element_text(size=12, face="bold"), legend.position = "right")
    print(p)
  }
  dev.off()
  
  # --- D. Output 2: Single PNG Grid ---
  png_file <- paste0("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_bars_", bin_name, ".png")
  
  p_grid <- ggplot(binned_data, aes(x = Bin_Label, y = Bin_Coverage_Mb, fill = Order)) +
    geom_bar(stat = "identity", width = (CURRENT_BIN_SIZE/1000000) * 0.9) +
    facet_wrap(~Chr_Label, ncol = 4, scales = "free_x") +
    scale_fill_manual(values = te_colors) +
    scale_y_continuous(limits = c(0, GLOBAL_MAX_Y), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    
    labs(title = paste0("Genomic Landscape of Transposable Elements (", bin_name, " resolution)"),
         subtitle = paste0("Stacked Coverage per ", bin_name, " bin"),
         x = "Position (Mb)", y = "Coverage (Mb)") +
    
    theme_bw() +
    theme(panel.grid = element_blank(), 
          strip.text = element_text(size=11, face="bold"), 
          strip.background = element_rect(fill="grey95"), 
          legend.position = "right")
  
  ggsave(png_file, plot = p_grid, width = 20, height = 15, dpi = 300)
}

print("All bin widths processed.")