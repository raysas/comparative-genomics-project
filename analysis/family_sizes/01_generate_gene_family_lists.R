library(dplyr)
library(readr)
library(stringr)

# Folder input
input_dir <- "../../output/clusters/"

target_files <- list(
  low  = "protein_families_filtered_blast_results_id30_qcov50_scov50_evalue1e-10_wcol12_network.tsv",
  high = "protein_families_filtered_blast_results_id50_qcov70_scov70_evalue1e-10_wcol12_network.tsv"
)

tsv_files <- list()

for (key in names(target_files)) {
  file_path <- file.path(input_dir, target_files[[key]])
  
  if (!file.exists(file_path)) {
    stop(paste("Required file not found:", file_path))
  }
  
  tsv_files[[key]] <- file_path
}

# Folder output
output_dir <- "../../output/family_sizes/genelists"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

for (level in names(tsv_files)) {
  input_file <- tsv_files[[level]]
  cat("\nProcessing", toupper(level), "stringency file:", basename(input_file), "\n")
  
  # Read TSV
  df <- read_tsv(input_file, show_col_types = FALSE)
  
  # Check required columns
  required_cols <- c("geneName", "family")
  if (!all(required_cols %in% colnames(df))) {
    stop(paste(
      "File", input_file, "missing required columns:", 
      paste(required_cols, collapse=", "),
      "Found:", paste(colnames(df), collapse=", ")
    ))
  }
  
  # Count family sizes
  family_sizes <- df %>%
    group_by(family) %>%
    summarise(size = n(), .groups = "drop")
  
  df <- df %>% left_join(family_sizes, by = "family")
  
  # Filter small (2–5) and large (>=6)
  small_genes <- df %>%
    filter(size >= 2, size <= 5) %>%
    arrange(geneName) %>%
    pull(geneName)
  
  large_genes <- df %>%
    filter(size >= 6) %>%
    arrange(geneName) %>%
    pull(geneName)
  
  # Output paths
  small_out <- file.path(output_dir, paste0("small_families_", level, "_stringency.txt"))
  large_out <- file.path(output_dir, paste0("large_families_", level, "_stringency.txt"))
  
  # Save outputs
  write_lines(small_genes, small_out)
  write_lines(large_genes, large_out)
  
  cat("Saved:", basename(small_out), ",", basename(large_out), "\n")
}

cat("\nALL DONE.\n")
