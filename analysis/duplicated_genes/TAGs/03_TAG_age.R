#!/usr/bin/env Rscript
# TAG age analysis - Compare Ks distributions
# Author: CG Project
# Date: 2025-11-23

library(tidyverse)

# Configuration
KS_FILE <- "../../output/statistics/ks_filtered.tsv"
TAG_FILE <- "../../output/statistics/TAGs_distance_100000bp.tsv"
GENE_ANNOTATION <- "../../output/statistics/gene_TAG_annotation.tsv"
DUPLICATION_TYPES <- "../../output/statistics/duplication_types.tsv"  # If available
OUTPUT_DIR <- "../../output/figures/"
STATS_DIR <- "../../output/statistics/"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STATS_DIR, recursive = TRUE, showWarnings = FALSE)

# === Load Data ===
cat("Loading data...\n")
ks_data <- read_tsv(KS_FILE, col_types = cols())
tag_pairs <- read_tsv(TAG_FILE, col_types = cols())
gene_annotation <- read_tsv(GENE_ANNOTATION, col_types = cols())

cat(sprintf("  Ks data: %d pairs\n", nrow(ks_data)))
cat(sprintf("  TAG pairs: %d pairs\n", nrow(tag_pairs)))
cat(sprintf("  Gene annotations: %d genes (%d TAGs)\n", 
            nrow(gene_annotation), sum(gene_annotation$Is_TAG)))

# === Merge TAG information with Ks ===
cat("\nMerging TAG annotation with Ks data...\n")

# Create TAG pair identifiers (need to match with Ks data)
# Assuming Ks data has Gene1 and Gene2 columns
tag_pair_ids <- tag_pairs %>%
  mutate(Pair_ID = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_"))

# Assuming ks_data has similar structure - adjust column names as needed
ks_with_tags <- ks_data %>%
  mutate(Pair_ID = paste(pmin(Gene1, Gene2), pmax(Gene1, Gene2), sep = "_"),
         Is_TAG = Pair_ID %in% tag_pair_ids$Pair_ID,
         Duplication_Type = ifelse(Is_TAG, "TAG", "Non-TAG"))

cat(sprintf("  TAG pairs with Ks: %d\n", sum(ks_with_tags$Is_TAG)))
cat(sprintf("  Non-TAG pairs: %d\n", sum(!ks_with_tags$Is_TAG)))

# === Statistical Comparison ===
cat("\n=== Statistical Tests: TAG vs Non-TAG Ks ===\n")

ks_tag <- ks_with_tags %>% filter(Is_TAG) %>% pull(Ks)
ks_nontag <- ks_with_tags %>% filter(!Is_TAG) %>% pull(Ks)

# Summary statistics
cat("TAG pairs:\n")
cat(sprintf("  N = %d\n", length(ks_tag)))
cat(sprintf("  Mean Ks = %.4f\n", mean(ks_tag)))
cat(sprintf("  Median Ks = %.4f\n", median(ks_tag)))
cat(sprintf("  SD = %.4f\n", sd(ks_tag)))

cat("\nNon-TAG pairs:\n")
cat(sprintf("  N = %d\n", length(ks_nontag)))
cat(sprintf("  Mean Ks = %.4f\n", mean(ks_nontag)))
cat(sprintf("  Median Ks = %.4f\n", median(ks_nontag)))
cat(sprintf("  SD = %.4f\n", sd(ks_nontag)))

# Mann-Whitney U test (non-parametric)
cat("\nMann-Whitney U test:\n")
wilcox_test <- wilcox.test(ks_tag, ks_nontag, alternative = "less")
print(wilcox_test)

cat(sprintf("\nConclusion: TAG Ks is %s than Non-TAG Ks (p = %.2e)\n",
            ifelse(wilcox_test$p.value < 0.05, "significantly lower", "not significantly different"),
            wilcox_test$p.value))

# T-test (parametric, for comparison)
cat("\nT-test (for comparison):\n")
t_test <- t.test(ks_tag, ks_nontag, alternative = "less")
print(t_test)

# Effect size (Cohen's d)
pooled_sd <- sqrt(((length(ks_tag)-1)*var(ks_tag) + (length(ks_nontag)-1)*var(ks_nontag)) / 
                  (length(ks_tag) + length(ks_nontag) - 2))
cohens_d <- (mean(ks_tag) - mean(ks_nontag)) / pooled_sd
cat(sprintf("\nEffect size (Cohen's d): %.3f\n", cohens_d))

# === Visualization ===
cat("\nGenerating plots...\n")

# Boxplot comparison
p1 <- ggplot(ks_with_tags, aes(x = Duplication_Type, y = Ks, fill = Duplication_Type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_violin(alpha = 0.3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("TAG" = "#E74C3C", "Non-TAG" = "#3498DB")) +
  labs(title = "Ks Distribution: TAGs vs Non-TAGs",
       subtitle = sprintf("Mann-Whitney U: p = %.2e | Cohen's d = %.3f", 
                         wilcox_test$p.value, cohens_d),
       x = "Duplication Type",
       y = "Ks (synonymous substitutions per site)") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold")) +
  annotate("text", x = 1, y = max(ks_with_tags$Ks)*0.9, 
           label = sprintf("n = %d\nMedian = %.3f", length(ks_tag), median(ks_tag))) +
  annotate("text", x = 2, y = max(ks_with_tags$Ks)*0.9, 
           label = sprintf("n = %d\nMedian = %.3f", length(ks_nontag), median(ks_nontag)))

ggsave(file.path(OUTPUT_DIR, "ks_TAG_vs_nonTAG_boxplot.pdf"), p1, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "ks_TAG_vs_nonTAG_boxplot.png"), p1, width = 8, height = 6, dpi = 300)

# Density plot comparison
p2 <- ggplot(ks_with_tags, aes(x = Ks, fill = Duplication_Type)) +
  geom_density(alpha = 0.6) +
  geom_vline(data = ks_with_tags %>% group_by(Duplication_Type) %>% 
             summarize(median_ks = median(Ks)),
             aes(xintercept = median_ks, color = Duplication_Type), 
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("TAG" = "#E74C3C", "Non-TAG" = "#3498DB")) +
  scale_color_manual(values = c("TAG" = "#E74C3C", "Non-TAG" = "#3498DB")) +
  labs(title = "Ks Distribution Comparison",
       subtitle = "TAGs show younger duplication ages",
       x = "Ks (synonymous substitutions per site)",
       y = "Density",
       fill = "Type",
       color = "Median") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "ks_TAG_vs_nonTAG_density.pdf"), p2, width = 10, height = 6)
ggsave(file.path(OUTPUT_DIR, "ks_TAG_vs_nonTAG_density.png"), p2, width = 10, height = 6, dpi = 300)

# === Orientation Analysis (if strand info available) ===
if ("Same_Orientation" %in% colnames(tag_pairs)) {
  cat("\n=== TAG Orientation Analysis ===\n")
  
  orientation_counts <- table(tag_pairs$Same_Orientation)
  cat("Orientation distribution:\n")
  print(orientation_counts)
  cat(sprintf("Same orientation: %d / %d (%.1f%%)\n",
              orientation_counts["TRUE"], sum(orientation_counts),
              orientation_counts["TRUE"]/sum(orientation_counts)*100))
  
  # Chi-square test for orientation bias
  # H0: orientation is random (50% same, 50% opposite)
  chisq_test <- chisq.test(orientation_counts, p = c(0.5, 0.5))
  cat("\nChi-square test (vs random orientation):\n")
  print(chisq_test)
  
  # Plot orientation
  p3 <- ggplot(tag_pairs, aes(x = Same_Orientation, fill = Same_Orientation)) +
    geom_bar() +
    geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5) +
    scale_fill_manual(values = c("FALSE" = "#95A5A6", "TRUE" = "#27AE60"),
                     labels = c("Opposite", "Same")) +
    labs(title = "TAG Gene Orientation",
         subtitle = sprintf("χ² test: p = %.2e", chisq_test$p.value),
         x = "Orientation",
         y = "Number of TAG pairs",
         fill = "Strand") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(OUTPUT_DIR, "TAG_orientation.pdf"), p3, width = 8, height = 6)
  ggsave(file.path(OUTPUT_DIR, "TAG_orientation.png"), p3, width = 8, height = 6, dpi = 300)
}

# === Save Results ===
cat("\nSaving results...\n")

# Statistical test results
test_results <- data.frame(
  Test = c("Mann-Whitney U", "T-test", "Cohen's d"),
  Statistic = c(wilcox_test$statistic, t_test$statistic, cohens_d),
  P_value = c(wilcox_test$p.value, t_test$p.value, NA),
  Conclusion = c(
    ifelse(wilcox_test$p.value < 0.05, "TAGs younger", "No difference"),
    ifelse(t_test$p.value < 0.05, "TAGs younger", "No difference"),
    ifelse(abs(cohens_d) > 0.5, "Medium/Large effect", "Small effect")
  )
)

write_tsv(test_results, file.path(STATS_DIR, "TAG_age_tests.tsv"))

# Summary statistics
summary_stats <- ks_with_tags %>%
  group_by(Duplication_Type) %>%
  summarize(
    N = n(),
    Mean_Ks = mean(Ks),
    Median_Ks = median(Ks),
    SD_Ks = sd(Ks),
    Q25 = quantile(Ks, 0.25),
    Q75 = quantile(Ks, 0.75)
  )

write_tsv(summary_stats, file.path(STATS_DIR, "TAG_ks_summary.tsv"))

cat("\n=== Analysis Complete ===\n")
cat(sprintf("Plots saved to: %s\n", OUTPUT_DIR))
cat(sprintf("Statistics saved to: %s\n", STATS_DIR))
