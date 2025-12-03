library(dplyr)
library(readr)

ks_table <- read_delim(
  "../../output/ks_results_low/gene_level_classification_with_ks_w_low.tsv",
  delim = "\t",
  col_types = cols()
)

mart <- read_delim(
  "../../data/mart_export.txt",
  delim = "\t",
  col_types = cols()
)

# Rename columns for clarity
mart <- mart %>%
  rename(
    GLYMA_ID = `Gene stable ID`,
    transcript = `Transcript stable ID`,
    GO_accession = `GO term accession`,
    GO_name = `GO term name`
  )

# Natural join on transcript ID
merged <- ks_table %>%
  left_join(mart, by = c("gene" = "transcript"))


merged <- merged %>%
  select(gene, GLYMA_ID, everything())


write.table(
  merged,
  "../../output/ks_results_low/ks_with_GLYMA_ID_GO_low.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

