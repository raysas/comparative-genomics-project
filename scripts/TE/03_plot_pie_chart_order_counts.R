# ------------------------------------------------------------------
# Script Name: 03_plot_pie_chart_order_counts.R
# Description: Generates a Pie Chart of TE Orders in Glycine max.
#              Class I elements are colored in Purple gradients.
#              Class II elements are colored in Blue gradients.
#              Data is sorted by abundance.
# ------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scales) 

# 1. Load Data 
# ----------------
data <- read.table("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_classes_and_orders_counts.txt", header=FALSE)
colnames(data) <- c("Count", "Class_Group", "Order_Name")

# 2. Process Data 
# ----------------
# Calculate Total for the Title
total_count <- sum(data$Count)

# Calculate fractions and percentages
data <- data %>%
  mutate(fraction = Count / sum(Count),
         percentage = round(fraction * 100, 1)) %>%
  # Sort by Class, then by Count Descending
  arrange(Class_Group, desc(Count)) 

# Create the Label Text: "Count (Percentage%)"
data$Label_Text <- paste0(comma(data$Count), " (", data$percentage, "%)")

# 3. Dynamic Color Generation
# ----------------
get_gradient <- function(colors_vec, n) {
  if (n == 0) return(c())
  colorRampPalette(colors_vec)(n)
}

# Extract names (Sorted by abundance)
class1_orders <- data$Order_Name[data$Class_Group == "Class_I"]
class2_orders <- data$Order_Name[data$Class_Group == "Class_II"]
unk_orders    <- data$Order_Name[data$Class_Group == "Unknown"]

# Assign Class I Colors (Warm Purples)
cols_c1 <- get_gradient(c("#7A0177", "#FA9FB5"), length(class1_orders))
names(cols_c1) <- class1_orders

# Assign Class II Colors (Cool Blues)
cols_c2 <- get_gradient(c("#084594", "#9ECAE1"), length(class2_orders))
names(cols_c2) <- class2_orders

# Assign Unknown Colors (Grey)
cols_unk <- rep("#999999", length(unk_orders))
names(cols_unk) <- unk_orders

# Combine into one master palette
full_palette <- c(cols_c1, cols_c2, cols_unk)

# 4. Plotting
# ----------------
# Lock factor levels to ensure correct sorting in the legend
data$Order_Name <- factor(data$Order_Name, levels = c(class1_orders, class2_orders, unk_orders))

p <- ggplot(data, aes(x="", y=Count, fill=Order_Name)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=full_palette) +
  
  # Label Logic: Show text only if slice is > 1.5% to prevent overlap
  geom_text(aes(label = ifelse(percentage > 1.5, Label_Text, "")), 
            position = position_stack(vjust = 0.5), 
            size=3, color="white", fontface="bold") +
  
  theme_void() +
  
  labs(title=paste0("TE Distribution by Order (Total n = ", comma(total_count), ")"), 
       subtitle="Class I: Retrotransposons (Purples) vs. Class II: DNA Transposons (Blues)",
       fill="TE Order")

# Save the plot
ggsave("C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_order_pie_chart.png", plot=p, width=8, height=7, dpi=300, bg = "white")