# ------------------------------------------------------------------
# Script Name: 10_plot_coverage_and_nesting_stats.R
#
# Description: 
#   PART 2: VISUALIZATION
#   Plots the statistics extracted in Script 09.
#   - PLOT 1: True Genomic Coverage (% of Genome) by Superfamily.
#   - PLOT 2: Nesting Index (Simple Sum / True BP) by Superfamily.
#     NOTE: This index specifically measures "Self-Nesting": the degree to 
#     which members of the same superfamily are stacked inside each other.
#   - LEGEND: Ordered by Class and Descending Abundance.
#
# Input:  /output/TE_true_coverage_stats.csv
# Output: /output/TE_true_coverage_bar.png
#         /output/TE_nesting_index_dot.png
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)

# ==================================================================
# 1. LOAD DATA
# ==================================================================
input_csv <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_true_coverage_stats.csv"
stats_data <- read.csv(input_csv, stringsAsFactors = FALSE)

# Extract Global Totals for Subtitles
total_stats <- stats_data %>% filter(Category == "Total")
total_pct <- total_stats$Genome_Pct
total_mb  <- total_stats$Mb

# Filter for Superfamily level data only
plot_data <- stats_data %>%
  filter(Category == "Superfamily") %>%
  # Filter out tiny families (< 0.01% genome) for cleaner plots
  filter(Genome_Pct > 0.01)

# ==================================================================
# 2. SETUP COLORS & ORDERING
# ==================================================================
# A. Sort Superfamilies (Bars)
plot_data <- plot_data %>% arrange(Class_Group, desc(Genome_Pct))

c1_sfs <- plot_data %>% filter(Class_Group == "Class I (Retro)") %>% pull(Type)
c2_sfs <- plot_data %>% filter(Class_Group == "Class II (DNA)") %>% pull(Type)
unk_sfs <- plot_data %>% filter(Class_Group == "Unknown") %>% pull(Type)
plot_data$Type <- factor(plot_data$Type, levels = c(c1_sfs, c2_sfs, unk_sfs))

# B. Sort Orders (Legend) - Descending Coverage
order_sums <- plot_data %>%
  group_by(Class_Group, Order) %>%
  summarise(Total_Cov = sum(Genome_Pct), .groups = "drop") %>%
  arrange(Class_Group, desc(Total_Cov))

c1_orders <- order_sums %>% filter(Class_Group == "Class I (Retro)") %>% pull(Order)
c2_orders <- order_sums %>% filter(Class_Group == "Class II (DNA)") %>% pull(Order)
unk_orders <- order_sums %>% filter(Class_Group == "Unknown") %>% pull(Order)

final_legend_levels <- c(c1_orders, c2_orders, unk_orders)
plot_data$Order <- factor(plot_data$Order, levels = final_legend_levels)

# C. Color Palette
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
# 3. PLOT A: TRUE COVERAGE (Biomass)
# ==================================================================
p_cov <- ggplot(plot_data, aes(x = Type, y = Genome_Pct, fill = Order)) +
  geom_bar(stat = "identity") +
  
  geom_text(aes(label = paste0(round(Genome_Pct, 2), "%")), 
            hjust = -0.1, size = 3, fontface="bold") +
  
  scale_fill_manual(values = te_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  
  coord_flip() +
  facet_wrap(~Class_Group, scales = "free_y", ncol = 1) + 
  
  theme_bw() +
  labs(title = "True Genomic Coverage by TE Superfamily",
       subtitle = paste0("Total Physical Coverage: ", total_pct, "% (", comma(total_mb), " Mb)"),
       x = "Superfamily",
       y = "% of Genome Covered (Physical Space)") +
  
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 8),
        legend.position = "right",
        panel.grid.major.y = element_blank())

ggsave("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_true_coverage_bar.png", 
       plot = p_cov, width = 12, height = 12, dpi = 300)

# ==================================================================
# 4. PLOT B: NESTING RATIOS
# ==================================================================
nest_data <- plot_data %>% filter(Genome_Pct > 0.1)

p_nest <- ggplot(nest_data, aes(x = Nesting_Index, y = Type, color = Order)) +
  
  geom_segment(aes(x = 1, xend = Nesting_Index, y = Type, yend = Type), color = "grey") +
  
  geom_point(size = 4) +
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  
  facet_wrap(~Class_Group, scales = "free_y", ncol = 1) +
  scale_color_manual(values = te_colors) +
  
  theme_bw() +
  labs(title = "Internal Nesting Analysis (Self-Stacking Index)",
     subtitle = "A ratio > 1.0 indicates 'Self-Nesting' where elements jump into members of their own superfamily.",
     x = "Nesting Index (Simple Sum / True Sum)",
     y = "Superfamily") +
  
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 9, face="bold"),
        legend.position = "right")

ggsave("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_nesting_index_dot.png", 
       plot = p_nest, width = 10, height = 10, dpi = 300)

print("Plots saved: TE_true_coverage_bar.png and TE_nesting_index_dot.png")