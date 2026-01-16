# ------------------------------------------------------------------
# Script Name: 06_plot_te_lengths_distribution.R
#
# Description: 
#   Visualizes the Length Distribution of TEs to assess genome integrity.
#   - USES: Density Plots (Smooth Histograms) on a Log10 Scale.
#   - LABELS: Labels peaks AND includes a sorted legend.
#   - COLORS: Matches BIOMASS-based gradients (Dark = Largest Total Size).
#   - SCALES: Fixed Y-axis to allow direct comparison between panels.
#
# Input:  /data/TEAnnotationFinal.gff3
# Output: /output/TE_length_distribution.png
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales)

# ==================================================================
# 1. DATA LOADING & PRE-PROCESSING
# ==================================================================
input_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/TEAnnotationFinal.gff3"

gff_data <- read.table(input_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)
raw_data <- gff_data %>% select(V3, V9)
colnames(raw_data) <- c("Full_Name", "Attributes")

# ==================================================================
# 2. PARSE HIERARCHY & EXTRACT LENGTHS
# ==================================================================
plot_data_prep <- raw_data %>%
  mutate(
    Length_Str = str_extract(Attributes, "Length=[0-9]+"), 
    Length = as.numeric(gsub("Length=", "", Length_Str))
  ) %>%
  separate(Full_Name, into = c("Class_Raw", "Order", "Superfamily"), sep = "/", fill = "right") %>%
  mutate(
    Class_Panel = case_when(
      grepl("Class II", Class_Raw) ~ "Class II (DNA Transposons)",
      grepl("Class I", Class_Raw) ~ "Class I (Retrotransposons)",
      TRUE ~ "Unknown"
    ),
    # Rename "-" or missing orders to "Unknown"
    Order = ifelse(Order == "-" | is.na(Order) | Order == "", "Unknown", Order)
  )

# Filter NA lengths and tiny fragments
plot_data <- plot_data_prep %>%
  filter(!is.na(Length))

# ==================================================================
# 3. SORTING & COLOR GENERATION (BY BIOMASS)
# ==================================================================
# Sort by TOTAL SIZE (Coverage) to match the Pie Chart colors.

order_size <- plot_data %>%
  group_by(Class_Panel, Order) %>%
  summarise(Total_BP = sum(Length), .groups = 'drop') %>%
  arrange(Class_Panel, desc(Total_BP)) 

# Extract sorted lists
c1_orders_sorted <- order_size %>% 
  filter(Class_Panel == "Class I (Retrotransposons)") %>% 
  pull(Order)

c2_orders_sorted <- order_size %>% 
  filter(Class_Panel == "Class II (DNA Transposons)") %>% 
  pull(Order)

unk_orders_sorted <- order_size %>% 
  filter(Class_Panel == "Unknown") %>% 
  pull(Order)

# Define Master Factor Levels
final_levels <- c(c1_orders_sorted, c2_orders_sorted, unk_orders_sorted)
plot_data$Order <- factor(plot_data$Order, levels = final_levels)

# Generate Gradients
get_gradient <- function(colors_vec, n) {
  if (n == 0) return(c())
  colorRampPalette(colors_vec)(n)
}

# Purples for Class I
cols_c1 <- get_gradient(c("#7A0177", "#FA9FB5"), length(c1_orders_sorted))
if(length(cols_c1) > 0) names(cols_c1) <- c1_orders_sorted

# Blues for Class II
cols_c2 <- get_gradient(c("#084594", "#9ECAE1"), length(c2_orders_sorted))
if(length(cols_c2) > 0) names(cols_c2) <- c2_orders_sorted

# Grey for Unknown
cols_unk <- rep("#999999", length(unk_orders_sorted))
if(length(cols_unk) > 0) names(cols_unk) <- unk_orders_sorted

te_colors <- c(cols_c1, cols_c2, cols_unk)

# ==================================================================
# 4. CALCULATE PEAKS FOR LABELING
# ==================================================================
label_coords <- plot_data %>%
  group_by(Class_Panel, Order) %>%
  filter(n() > 1) %>%
  summarise(
    Peak_Log_X = density(log10(Length))$x[which.max(density(log10(Length))$y)],
    Peak_Y = max(density(log10(Length))$y),
    .groups = 'drop'
  ) %>%
  mutate(Peak_X = 10^Peak_Log_X)

# ==================================================================
# 5. VISUALIZATION
# ==================================================================
p <- ggplot(plot_data, aes(x = Length, color = Order, fill = Order)) +
  
  geom_density(alpha = 0.2, size = 0.8) +
  
  # Peak Labels
  geom_label(data = label_coords, 
             aes(x = Peak_X, y = Peak_Y, label = Order),
             color = "black", fill = "white", alpha = 0.8, size = 3,
             fontface = "bold",
             show.legend = FALSE) +
  
  scale_x_log10(breaks = c(10, 50, 100, 500, 1000, 5000, 10000, 15000),
                labels = c("10bp", "50bp", "100bp", "500bp", "1kb", "5kb", "10kb", "15kb")) +
  
  # Scales = "fixed" ensures Y-axis is identical across panels
  facet_wrap(~Class_Panel, scales = "fixed", ncol = 1) +
  
  scale_color_manual(values = te_colors) +
  scale_fill_manual(values = te_colors) +
  
  theme_bw() +
  
  labs(title = "TE Length Distribution Profile",
       subtitle = "Density plot showing size frequency (Log-scaled X-axis)",
       x = "Length (bp)",
       y = "Density (Frequency)",
       fill = "Order",
       color = "Order") + 
  
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.position = "right", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# ==================================================================
# 6. SAVE OUTPUT
# ==================================================================
output_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_length_distribution.png"
ggsave(output_path, plot = p, width = 12, height = 10, dpi = 300)

print(paste("Plot saved to:", output_path))