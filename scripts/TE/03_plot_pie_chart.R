# ------------------------------------------------------------------
# Script Name: 03_plot_pie_chart_sorted.R
# Description: Generates a Pie Chart of TE Orders.
#              Class I = Warm Purples, Class II = Cool Blues.
#              SORTED by Abundance (Darkest Color = Largest Slice).
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 1. Load Data 
# ----------------
data <- read.table("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_classes_and_orders_counts.txt", header=FALSE)
colnames(data) <- c("Count", "Class_Group", "Order_Name")

# 2. Process Data (THIS IS THE ONLY MAJOR CHANGE)
# -----------------------------------------------
# Calculate fractions and percentages for labels
data <- data %>%
  mutate(fraction = Count / sum(Count),
         percentage = round(fraction * 100, 1)) %>%
  # CHANGE: Sort by Class, then by Count Descending (Biggest First)
  # This ensures the biggest groups get the darkest colors
  arrange(Class_Group, desc(Count)) 

# 3. Dynamic Color Generation
# ----------------
get_gradient <- function(colors_vec, n) {
  if (n == 0) return(c())
  colorRampPalette(colors_vec)(n)
}

# -- Extract names (These are now sorted Big -> Small) --
class1_orders <- data$Order_Name[data$Class_Group == "Class_I"]
class2_orders <- data$Order_Name[data$Class_Group == "Class_II"]
unk_orders    <- data$Order_Name[data$Class_Group == "Unknown"]

# -- Assign Class I Colors (Warm Purples) --
# First item (Biggest) gets #7A0177 (Dark Purple)
cols_c1 <- get_gradient(c("#7A0177", "#FA9FB5"), length(class1_orders))
names(cols_c1) <- class1_orders

# -- Assign Class II Colors (Cool Blues) --
# First item (Biggest) gets #084594 (Dark Blue)
cols_c2 <- get_gradient(c("#084594", "#9ECAE1"), length(class2_orders))
names(cols_c2) <- class2_orders

# -- Assign Unknown Colors (Grey) --
cols_unk <- rep("#999999", length(unk_orders))
names(cols_unk) <- unk_orders

# Combine into one master palette
full_palette <- c(cols_c1, cols_c2, cols_unk)

# 4. Plotting
# ----------------
# Lock the factor levels to the sorted order so ggplot draws them correctly
data$Order_Name <- factor(data$Order_Name, levels = c(class1_orders, class2_orders, unk_orders))

p <- ggplot(data, aes(x="", y=Count, fill=Order_Name)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=full_palette) +
  
  # Label Logic: Only show text if the slice is bigger than 1%
  geom_text(aes(label = ifelse(percentage > 1, paste0(percentage, "%"), "")), 
            position = position_stack(vjust = 0.5), size=3.5, color="white", fontface="bold") +
  
  theme_void() +
  labs(title="TE Distribution by Order", 
       subtitle="Class I (Purples) vs. Class II (Blues)",
       fill="TE Order")

# Save the plot with white background
ggsave("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_order_pie_chart.png", plot=p, width=7, height=6, dpi=300, bg = "white")