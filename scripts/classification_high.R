library(dplyr)
library(tidyr)


# 1. LOAD PAIRWISE DATA

data <- read.table(
  "../output/ks_results_high/ks_results_high.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

data <- data %>% select(gene1, gene2, ks, ka_ks)


# 2. LONG FORMAT (one row per gene)

long <- data %>%
  pivot_longer(cols = c(gene1, gene2),
               names_to = "pos",
               values_to = "gene") %>%
  select(gene, ks, ka_ks)



# 3. AGE CLASSIFICATION (using ks_min)
#    ks_min < 2  → Recent
#    ks_min >= 2 → Old


age_info <- long %>%
  group_by(gene) %>%
  summarise(
    ks_used = min(ks, na.rm = TRUE),
    age_class = ifelse(ks_used >= 2, "Old", "Recent")
  )


# 4. SELECTION CLASSIFICATION (using w_min / w_max)
#    If w_min > 1 → Positive
#    If w_max < 1 → Purifying
#    Else         → Neutral


sel_info <- long %>%
  group_by(gene) %>%
  summarise(
    w_min = min(ka_ks, na.rm = TRUE),
    w_max = max(ka_ks, na.rm = TRUE)
  ) %>%
  mutate(
    selection_class = case_when(
      w_min > 1 ~ "Positive_selection",
      w_max < 1 ~ "Purifying_selection",
      TRUE      ~ "Neutral"
    ),
    # KEEP ONLY the determining value
    w_used = case_when(
      selection_class == "Positive_selection" ~ w_min,
      selection_class == "Purifying_selection" ~ w_max,
      selection_class == "Neutral" ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  select(gene, selection_class, w_used)



# 5. MERGE AGE + SELECTION INTO FINAL TABLE

final_gene_table <- age_info %>%
  left_join(sel_info, by = "gene") %>%
  select(gene, age_class, ks_used, selection_class, w_used)


# 6. SAVE OUTPUT

out_file <- "../output/ks_results_high/gene_level_classification_with_ks_w_high.tsv"

write.table(final_gene_table, out_file,
            sep="\t", quote=FALSE, row.names=FALSE)

cat(" Final gene-level table saved to:", out_file, "\n")
