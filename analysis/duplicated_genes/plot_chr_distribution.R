library(ggplot2)
library(ggpattern)
library(dplyr)
library(tidyr)

dup_df<-read.csv('../../output/info/duplicated_genes_info_id50_qcov70_scov70_evalue1e-10_wcol12.csv')
prot_df<-read.csv('../../output/info/protein_info_longest.csv')

total_counts <- prot_df %>% count(chromosome, name = "total_genes")
dup_counts <- dup_df %>% count(chromosome, name = "duplicated_genes")

# ------------ dup and total genes per chromosome -----------------

merged_df <- total_counts %>%
  left_join(dup_counts, by = "chromosome") %>%
  mutate(
    duplicated_genes = ifelse(is.na(duplicated_genes), 0, duplicated_genes),
    non_duplicated = total_genes - duplicated_genes
  )

plot_df <- merged_df %>%
  pivot_longer(cols = c(duplicated_genes, non_duplicated),
               names_to = "gene_type",
               values_to = "count")

ggplot(plot_df, aes(x = chromosome, y = count, fill = chromosome, pattern = gene_type)) +
  geom_bar_pattern(
    stat = "identity",
    width = 0.8,
    color = "white",
    linewidth = 0.3,
    pattern_fill = "black",
    pattern_color = "black",
    pattern_density = 0.15,
    pattern_spacing = 0.03,
    pattern_key_scale_factor = 0.5
  ) +
  scale_fill_manual(values = chrom_colors) +
  scale_pattern_manual(values = c("duplicated_genes" = "none", "non_duplicated" = "stripe"),
                       labels = c("duplicated genes", "non-duplicated")) +
  labs(title = "total genes vs duplicated genes per chromosome",
       x = "chromosome", y = "gene count",
       fill = "chromosome", pattern = "gene type") +
  base_theme +
  custom_theme

# -------------- tag and non tag dup genes per chromosome -----------------



families_tag_dup_df <- dup_df %>%
  inner_join(families_df, by = c("peptide_id"="geneName")) %>%
  left_join(tags_df %>% select(peptide_id, TAG), by = "peptide_id") %>%
  mutate(
    is_TAG = ifelse(!is.na(TAG) & TAG != 0, TRUE, FALSE),
    TAG_id = ifelse(is_TAG, paste0(family, "_TAG", TAG), NA)
  )

# -- plot prep
tag_dup_counts <- families_tag_dup_df %>%
  group_by(chromosome, is_TAG) %>%
  summarise(count = n()) %>%
  ungroup()
tag_dup_counts$is_TAG <- ifelse(tag_dup_counts$is_TAG, "TAG duplicated genes", "non-TAG duplicated genes")

ggplot(tag_dup_counts, aes(x = chromosome, y = count, fill = is_TAG)) +
  geom_bar(stat = "identity",
           width = 0.8,
           color = "white",
           linewidth = 0.3) +
  scale_fill_manual(values = chrom_colors) +
  labs(title = "TAG vs non-TAG duplicated genes per chromosome",
       x = "chromosome", y = "gene count",
       fill = "duplicated gene type") +
  base_theme +
  custom_theme

