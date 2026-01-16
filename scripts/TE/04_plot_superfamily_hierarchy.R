# ------------------------------------------------------------------
# Script Name: 04_plot_superfamily_hierarchy.R
# Description: Generates a faceted bar chart visualizing Transposable 
#              Element (TE) superfamily abundance in Glycine max.
#
# Key Features:
#   - Facets data into two panels: Class I (Retrotransposons) and 
#     Class II (DNA Transposons).
#   - Groups superfamilies by their taxonomic Order.
#   - Sorts elements by abundance (descending).
#   - Enhances legend with specific Subclass information.
#   - Applies a distinct color gradient (Purples for Class I, Blues for Class II).
#
# Input:  /output/TE_classification_counts.txt (Tab-separated counts)
# Output: /output/TE_superfamily_hierarchy.png
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# 1. Data Loading and Cleaning
# ------------------------------------------------------------------
# Reads the raw text file containing TE counts and classification strings.
input_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_classification_counts.txt"

raw_lines <- readLines(input_path)
data <- data.frame(raw_text = trimws(raw_lines), stringsAsFactors = FALSE)

# Separates the count (numeric) from the classification string.
# 'extra = "merge"' ensures complex names with spaces are preserved.
clean_data <- data %>%
  separate(raw_text, into = c("Count", "Full_Name"), sep = "\\s+", extra = "merge") %>%
  mutate(Count = as.numeric(Count))

# 2. Hierarchy Parsing and Classification
# ------------------------------------------------------------------
# Splits the full classification string into Class, Order, and Superfamily.
# Defines panel groups and identifies Subclasses for legend labeling.
plot_data_prep <- clean_data %>%
  separate(Full_Name, into = c("Class_Raw", "Order", "Superfamily"), sep = "/", fill = "right") %>%
  mutate(
    # Define Panel Groups (Panel Headers)
    Class_Group = case_when(
      grepl("Class II", Class_Raw) ~ "Class II (DNA Transposons)",
      grepl("Class I", Class_Raw) ~ "Class I (Retrotransposons)",
      TRUE ~ "Unknown"
    ),
    # Standardize missing superfamily names
    Superfamily = ifelse(Superfamily == "-" | is.na(Superfamily) | Superfamily == "", "Unspecified", Superfamily)
  )

# 3. Filtering
# ------------------------------------------------------------------
# Removes unclassified entries. 
plot_data <- plot_data_prep %>%
  filter(Class_Group != "Unknown") %>%
  # Generates unique y-axis labels to prevent merging of "Unknown" categories across Orders
  mutate(
    Plot_Label = ifelse(Superfamily == "Unknown" | Superfamily == "Unspecified", 
                        paste0("Unknown (", Order, ")"), 
                        Superfamily)
  )

# 4. Abundance Calculation and Ordering
# ------------------------------------------------------------------
# Calculates total abundance per Order to determine the sort order.
order_stats <- plot_data %>%
  group_by(Class_Group, Order) %>%
  summarise(Total_Abundance = sum(Count), .groups = 'drop') 

# Extracts sorted lists of Orders for the Legend (Highest abundance first).
c1_legend <- order_stats %>%
  filter(Class_Group == "Class I (Retrotransposons)" & Order != "Unknown") %>%
  arrange(desc(Total_Abundance)) %>%
  pull(Order)

c2_legend <- order_stats %>%
  filter(Class_Group == "Class II (DNA Transposons)" & Order != "Unknown") %>%
  arrange(desc(Total_Abundance)) %>%
  pull(Order)

# Defines the master factor levels for the legend.
final_legend_levels <- c(c1_legend, c2_legend, "Unknown")
plot_data$Order <- factor(plot_data$Order, levels = final_legend_levels)

# 5. Label Generation (Subclass Detection)
# ------------------------------------------------------------------
# Creates a mapping vector to append "Subclass 1" or "Subclass 2" to legend text
# based on the raw classification string, without altering the underlying data structure.
label_map_df <- plot_data %>%
  select(Order, Class_Raw) %>%
  distinct() %>%
  group_by(Order) %>%
  filter(row_number() == 1) %>% 
  mutate(
    Suffix = case_when(
      grepl("subclass [1I]", Class_Raw, ignore.case = TRUE) ~ " (Subclass 1)",
      grepl("subclass 2", Class_Raw, ignore.case = TRUE) ~ " (Subclass 2)",
      TRUE ~ ""
    ),
    New_Label = paste0(Order, Suffix)
  )

legend_labels <- setNames(label_map_df$New_Label, label_map_df$Order)

# 6. Plot Data Sorting
# ------------------------------------------------------------------
# Sorts the data frame to ensure bars are stacked correctly in the plot.
plot_data_sorted <- plot_data %>%
  left_join(order_stats, by = c("Class_Group", "Order")) %>%
  mutate(Is_Unknown = ifelse(Order == "Unknown", 1, 0)) %>%
  arrange(
    Class_Group,          
    desc(Is_Unknown),     # Forces Unknowns to the bottom of the stack
    Total_Abundance,      # Groups by Order size
    Count                 # Sorts individual families
  )

# Locks the sort order for plotting
plot_data_sorted$Plot_Label <- factor(plot_data_sorted$Plot_Label, levels = unique(plot_data_sorted$Plot_Label))

# 7. Color Palette Generation
# ------------------------------------------------------------------
# Generates dynamic color gradients based on the number of Orders in each Class.
get_gradient <- function(colors_vec, n) {
  if (n == 0) return(c())
  colorRampPalette(colors_vec)(n)
}

# Class I: Purple Gradient
cols_c1 <- get_gradient(c("#7A0177", "#FA9FB5"), length(c1_legend))
if(length(cols_c1) > 0) names(cols_c1) <- c1_legend

# Class II: Blue Gradient
cols_c2 <- get_gradient(c("#084594", "#9ECAE1"), length(c2_legend))
if(length(cols_c2) > 0) names(cols_c2) <- c2_legend

cols_unk <- c("Unknown" = "#999999")
te_colors <- c(cols_c1, cols_c2, cols_unk)

# 8. Visualization (ggplot2)
# ------------------------------------------------------------------
# Calculate global max for consistent scaling across panels
global_max_count <- max(plot_data_sorted$Count)

p <- ggplot(plot_data_sorted, aes(x = Plot_Label, y = Count, fill = Order)) +
  
  geom_bar(stat = "identity") +
  
  # Revert to scales="free" to keep the nice separation of families (Vertical axis)
  # But use expand_limits to force the Count axis (Horizontal) to match
  facet_wrap(~Class_Group, scales = "free", ncol = 2) +
  
  coord_flip() +
  
  # Force both panels to span the full range of counts
  expand_limits(y = global_max_count) +
  
  # Apply colors and custom legend labels
  scale_fill_manual(values = te_colors, 
                    breaks = final_legend_levels,
                    labels = legend_labels[final_legend_levels]) +
  
  scale_y_continuous(labels = comma) +
  theme_bw() +
  
  labs(title = "TE Superfamily Abundance Distribution",
       subtitle = "Grouped by Order",
       x = "Superfamily",
       y = "Count",
       fill = "Order") +
  
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "right",
        # Remove grid lines for cleaner appearance
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# 9. Save Output
# ------------------------------------------------------------------
output_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_superfamily_hierarchy.png"
ggsave(output_path, plot = p, width = 14, height = 8, dpi = 300)

print(paste("Plot saved successfully to:", output_path))