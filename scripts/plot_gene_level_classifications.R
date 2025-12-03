library(dplyr)
library(ggplot2)
library(tidyr)




scientific_theme <- function() {
  theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text  = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10)
    )
}


# 1. LOAD FINAL GENE TABLE

df <- read.table(
  "../output/ks_results/gene_level_classification_with_ks_w.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)


# 2. CREATE OUTPUT DIR FOR FIGURES


plot_dir <- "../figures/ks"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# 3. PLOT: Histogram of Ks_used (gene ages)


p1 <- ggplot(df, aes(x = ks_used)) +
  geom_histogram(binwidth = 0.2, fill = "#4063D8", color = "white") +
  labs(title = "Distribution of synonymous substitution rate (Ks)",
       x = "Ks",
       y = "Gene count") +
  scientific_theme()

ggsave(file.path(plot_dir, "hist_ks_used.png"),
       p1, width = 7, height = 5, dpi = 600)



# 4. PLOT: Histogram of w_used (Ka/Ks)


p2 <- ggplot(df %>% filter(!is.na(w_used)), aes(x = w_used)) +
  geom_histogram(binwidth = 0.1, fill = "#CC79A7", color = "white") +
  labs(title = "Distribution of Selection Pressure (Ka/Ks)",
       x = "Ka/Ks",
       y = "Gene count") +
  scientific_theme()

ggsave(file.path(plot_dir, "hist_w_used.png"),
       p2, width = 7, height = 5, dpi = 600)



# 5. PLOT: Barplot of age_class


p3 <- ggplot(df, aes(x = age_class)) +
  geom_bar(fill = "#4DAF4A") +
  labs(title = "Gene Duplication Age Classification",
       x = "Duplication age",
       y = "Gene count") +
  scientific_theme()

ggsave(file.path(plot_dir, "bar_age_class.png"),
       p3, width = 6, height = 4, dpi = 600)



# 6. PLOT: Barplot of selection class


p4 <- ggplot(df, aes(x = selection_class)) +
  geom_bar(fill = "#984EA3") +
  labs(title = "Evolutionary Selection Class per Gene",
       x = "Selection class",
       y = "Gene count") +
  scientific_theme()

ggsave(file.path(plot_dir, "bar_selection_class.png"),
       p4, width = 7, height = 4, dpi = 600)



# 7. HEATMAP: Age class × Selection class


tab <- table(df$age_class, df$selection_class)
tab_df <- as.data.frame(tab)

p5 <- ggplot(tab_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "#e0ecf4", high = "#8856a7") +
  labs(title = "Relationship Between Duplication Age and Selection",
       x = "Selection class",
       y = "Duplication age") +
  scientific_theme()

ggsave(file.path(plot_dir, "heatmap_age_selection.png"),
       p5, width = 7, height = 5, dpi = 600)



# 8. SCATTERPLOT: Ks_used vs w_used

p6 <- ggplot(df %>% filter(!is.na(w_used)),
             aes(x = ks_used, y = w_used, color = selection_class)) +
  geom_point(alpha = 0.6, size = 1.7) +
  scale_color_manual(values = c(
    "Positive_selection"  = "#D55E00",
    "Purifying_selection" = "#009E73",
    "Neutral_or_mixed"    = "#888888"
  )) +
  labs(title = "Relationship Between Ks and Selection Pressure (Ka/Ks)",
       x = "Ks",
       y = "w = Ka/Ks (value used)") +
  scientific_theme()

ggsave(file.path(plot_dir, "scatter_ks_vs_w.png"),
       p6, width = 7, height = 5, dpi = 600)

cat("\n All plots saved in:", plot_dir, "\n")
