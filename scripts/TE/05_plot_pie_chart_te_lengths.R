# ------------------------------------------------------------------
# Script Name: 05_plot_pie_chart_te_lengths.R
#
# Description: 
#   Generates a Pie Chart representing the Genomic Coverage (Biomass) 
#   of Transposable Elements (TEs).
#
# Methodology:
#   1. Input Parsing: Reads GFF3 attributes to extract explicit TE lengths.
#   2. Aggregation: Sums the lengths (bp) for each TE Order.
#   3. Conversion: Converts raw base pairs to Megabases (Mb).
#   4. Visualization: Plots a Pie Chart with:
#      - Dynamic Colors (Purple for Retro, Blue for DNA).
#      - LABELS: Only shows numbers (Mb + %) for slices > 1%.
#      - LEGEND: Grouped by Class (Retro -> DNA -> Unknown) and sorted by size.
#
# Input:  /data/TEAnnotationFinal.gff3
# Output: /output/TE_coverage_pie_chart.png
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggrepel) # Required for non-overlapping labels

# ==================================================================
# 1. DATA LOADING & PRE-PROCESSING
# ==================================================================
input_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/TEAnnotationFinal.gff3"

gff_data <- read.table(input_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)
raw_data <- gff_data %>% select(V3, V9)
colnames(raw_data) <- c("Full_Name", "Attributes")

# ==================================================================
# 2. LENGTH EXTRACTION & HIERARCHY PARSING
# ==================================================================
plot_data_prep <- raw_data %>%
  mutate(
    # Extract Length from Attributes
    Length_Str = str_extract(Attributes, "Length=[0-9]+"), 
    Length = as.numeric(gsub("Length=", "", Length_Str))
  ) %>%
  separate(Full_Name, into = c("Class_Raw", "Order", "Superfamily"), sep = "/", fill = "right") %>%
  mutate(
    Class_Group = case_when(
      grepl("Class II", Class_Raw) ~ "Class II (DNA Transposons)",
      grepl("Class I", Class_Raw) ~ "Class I (Retrotransposons)",
      TRUE ~ "Unknown"
    ),
    Order = ifelse(Order == "-" | is.na(Order) | Order == "", "Unknown", Order)
  )

# ==================================================================
# 3. AGGREGATION (CALCULATE GENOME COVERAGE)
# ==================================================================
coverage_data <- plot_data_prep %>%
  filter(!is.na(Length)) %>%
  group_by(Class_Group, Order) %>%
  summarise(Total_BP = sum(Length), .groups = 'drop') %>%
  mutate(
    Total_Mb = round(Total_BP / 1000000, 2)
  ) %>%
  filter(Total_Mb > 0)

# Calculate Stats and Labels
total_te_space <- sum(coverage_data$Total_Mb)

coverage_final <- coverage_data %>%
  mutate(
    fraction = Total_BP / sum(Total_BP),
    percentage = round(fraction * 100, 1),
    
    # --- LABEL LOGIC: Only show numbers if > 1% ---
    Label_Text = ifelse(percentage > 1,
                        paste0(comma(Total_Mb), " Mb (", percentage, "%)"),
                        "") # Empty string for small slices
  ) %>%
  # Sort descending by size (Used for extracting ordered lists below)
  arrange(desc(Total_BP))

# ==================================================================
# 4. COLOR PALETTE & LEGEND ORDERING
# ==================================================================
# We need to extract the orders *sorted by size* within their classes.

# 1. Class I (Retro) - Sorted Big to Small
c1_orders <- coverage_final %>% 
  filter(Class_Group == "Class I (Retrotransposons)") %>% 
  pull(Order)

# 2. Class II (DNA) - Sorted Big to Small
c2_orders <- coverage_final %>% 
  filter(Class_Group == "Class II (DNA Transposons)") %>% 
  pull(Order)

# 3. Unknown
unk_orders <- coverage_final %>% 
  filter(Class_Group == "Unknown") %>% 
  pull(Order)

# Define the Master Legend Levels (Grouped by Class)
final_legend_levels <- c(c1_orders, c2_orders, unk_orders)

# GENERATE COLORS
get_gradient <- function(colors_vec, n) {
  if (n == 0) return(c())
  colorRampPalette(colors_vec)(n)
}

cols_c1 <- get_gradient(c("#7A0177", "#FA9FB5"), length(c1_orders))
if(length(cols_c1) > 0) names(cols_c1) <- c1_orders

cols_c2 <- get_gradient(c("#084594", "#9ECAE1"), length(c2_orders))
if(length(cols_c2) > 0) names(cols_c2) <- c2_orders

cols_unk <- rep("#999999", length(unk_orders))
if(length(cols_unk) > 0) names(cols_unk) <- unk_orders

full_palette <- c(cols_c1, cols_c2, cols_unk)

# ==================================================================
# 5. VISUALIZATION
# ==================================================================
# Apply the Grouped Factor Levels
coverage_final$Order <- factor(coverage_final$Order, levels = final_legend_levels)

p <- ggplot(coverage_final, aes(x = "", y = percentage, fill = Order)) +
  
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  
  scale_fill_manual(values = full_palette) +
  
  # Labels: Use ggrepel, will interpret "" as no label
  geom_text_repel(aes(label = Label_Text), 
                  position = position_stack(vjust = 0.5), 
                  size = 3, 
                  color = "white", 
                  fontface = "bold",
                  show.legend = FALSE,
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  
  theme_void() +
  
  labs(title = paste0("TE Total Lengths Distribution by Order (Total = ", comma(total_te_space), " Mb)"), 
       subtitle = "Class I: Retrotransposons (Purples) vs. Class II: DNA Transposons (Blues)",
       fill = "TE Order") +
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

# ==================================================================
# 6. SAVE OUTPUT
# ==================================================================
output_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_coverage_pie_chart.png"
ggsave(output_path, plot = p, width = 10, height = 10, dpi = 300, bg = "white")

print(paste("Plot saved successfully to:", output_path))