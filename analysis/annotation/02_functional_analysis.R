library(dplyr)
library(readr)
library(gprofiler2)
library(ggplot2)
library(forcats)

# ==========================================================
# 1. LOAD DATA (resolve project-root-relative paths)
# ==========================================================
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- NULL
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  fr <- sys.frames()
  if (!is.null(fr) && length(fr) > 0) {
    maybe <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
    if (!is.null(maybe)) script_path <- normalizePath(maybe)
  }
}

find_project_root <- function(start_dir) {
  dir <- normalizePath(start_dir)
  repeat {
    if (file.exists(file.path(dir, ".git")) ||
        file.exists(file.path(dir, "readme.md")) ||
        file.exists(file.path(dir, "README.md")) ||
        file.exists(file.path(dir, "scripts"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) break
    dir <- parent
  }
  return(NULL)
}

if (!is.null(script_path) && file.exists(script_path)) {
  script_dir <- dirname(script_path)
} else {
  script_dir <- normalizePath(".")
}
project_dir <- find_project_root(script_dir)
if (is.null(project_dir)) project_dir <- normalizePath(".")

in_file <- file.path(project_dir, "output", "ks_results", "ks_with_GLYMA_ID_GO.tsv")
if (!file.exists(in_file)) stop("Input file not found: ", in_file)

df <- read_delim(
  in_file,
  delim = "\t",
  col_types = cols()
)

out_enrich_dir <- file.path(project_dir, "output", "enrichment_gprofiler")
out_fig_dir <- file.path(project_dir, "figures", "enrichment_gprofiler")
dir.create(out_enrich_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

df <- df %>% 
  filter(!is.na(gene)) %>% 
  mutate(GLYMA_ID = as.character(GLYMA_ID))


# ==========================================================
# 2. GENE SETS
# ==========================================================
recent_genes  <- df %>% filter(age_class == "Recent") %>% pull(GLYMA_ID) %>% unique()
old_genes     <- df %>% filter(age_class == "Old")    %>% pull(GLYMA_ID) %>% unique()
purifying_genes <- df %>% filter(selection_class == "Purifying_selection") %>% pull(GLYMA_ID) %>% unique()
positive_genes  <- df %>% filter(selection_class == "Positive_selection") %>% pull(GLYMA_ID) %>% unique()

df <- df %>% mutate(ks_group = ifelse(ks_used <= 2, "<=2", ">2"))
ks_le2_genes <- df %>% filter(ks_group == "<=2") %>% pull(GLYMA_ID) %>% unique()
ks_gt2_genes <- df %>% filter(ks_group == ">2")  %>% pull(GLYMA_ID) %>% unique()


# ==========================================================
# 3. RUN gPROFILER
# ==========================================================
run_gprof <- function(genes, name){
  if (length(genes) == 0){
    message(paste("Skipping", name, ": empty list"))
    return(NULL)
  }
  
  gp <- gost(
    query = genes,
    organism = "gmax",
    sources = c("GO:BP", "GO:MF", "GO:CC")
  )
  
  if (!is.null(gp$result)){
    out_csv <- paste0("../../output/enrichment_gprofiler/", name, "_gprof.csv")

    # Robustly sanitize the result: convert any list-columns or multi-valued
    # elements into semicolon-separated strings and coerce all columns to
    # character so write.csv/write.table won't fail.
    sanitize_df <- function(x) {
      df <- as.data.frame(x, stringsAsFactors = FALSE)
      for (colname in names(df)) {
        col_vals <- df[[colname]]
        # detect list columns or elements with length > 1
        is_problem <- is.list(col_vals) || any(sapply(col_vals, function(z) is.list(z) || length(z) > 1))
        if (is_problem) {
          df[[colname]] <- sapply(col_vals, function(z) {
            if (is.null(z)) return(NA_character_)
            if (is.atomic(z)) return(paste(as.character(z), collapse = ";"))
            paste(unlist(z), collapse = ";")
          }, USE.NAMES = FALSE)
        } else {
          # ensure atomic columns are character
          df[[colname]] <- as.character(col_vals)
        }
      }
      df
    }

    res_df <- sanitize_df(gp$result)
    write.csv(res_df, out_csv, row.names = FALSE)
    return(res_df)
  }
  
  return(NULL)
}

# Run all analyses
res_recent    <- run_gprof(recent_genes,    "Recent")
res_old       <- run_gprof(old_genes,       "Old")
res_pur       <- run_gprof(purifying_genes, "Purifying")
res_pos       <- run_gprof(positive_genes,  "Positive")
res_kle2      <- run_gprof(ks_le2_genes,     "Ks_leq_2")
res_kgt2      <- run_gprof(ks_gt2_genes,     "Ks_gt_2")


# ==========================================================
# 4. VISUALISATION FUNCTIONS
# ==========================================================

# ---- A. Dotplot ----
plot_dot <- function(res, title, filename){
  df <- res %>% 
    filter(source == "GO:BP") %>% 
    # ensure p_value is numeric (g:Profiler results may be character after sanitization)
    mutate(p_value = as.numeric(p_value)) %>%
    filter(!is.na(p_value)) %>%
    arrange(p_value) %>% 
    slice_head(n = 10) %>% 
    mutate(term_name = fct_reorder(term_name, -log10(p_value)))
  
  p <- ggplot(df, aes(x = term_name, y = -log10(p_value))) +
    geom_point(size = 4, color = "#0072B2") +
    coord_flip() +
    theme_minimal(base_size = 13) +
    labs(
      title = title,
      x = "GO term (BP)",
      y = "-log10(p-value)"
    )
  
  ggsave(paste0("../../figures/enrichment_gprofiler/", filename, "_dotplot.png"), 
         p, width = 10, height = 6)
}

# ---- B. Barplot ----
plot_bar <- function(res, title, filename){
  df <- res %>% 
    filter(source == "GO:BP") %>% 
    mutate(p_value = as.numeric(p_value)) %>%
    filter(!is.na(p_value)) %>%
    arrange(p_value) %>% 
    slice_head(n = 10) %>% 
    mutate(term_name = fct_reorder(term_name, -log10(p_value)))
  
  p <- ggplot(df, aes(x = term_name, y = -log10(p_value))) +
    geom_bar(stat = "identity", fill = "#E69F00") +
    coord_flip() +
    theme_minimal(base_size = 13) +
    labs(
      title = title,
      x = "GO term (BP)",
      y = "-log10(p-value)"
    )
  
  ggsave(paste0("../../figures/enrichment_gprofiler/", filename, "_barplot.png"), 
         p, width = 10, height = 6)
}

# ==========================================================
# 5. GENERATE ALL PLOTS
# ==========================================================
if (!is.null(res_recent))    { plot_dot(res_recent, "Recent duplicates", "recent");       plot_bar(res_recent, "Recent duplicates", "recent") }
if (!is.null(res_old))       { plot_dot(res_old, "Old duplicates", "old");               plot_bar(res_old, "Old duplicates", "old") }
if (!is.null(res_pur))       { plot_dot(res_pur, "Purifying selection", "purifying");     plot_bar(res_pur, "Purifying selection", "purifying") }
if (!is.null(res_pos))       { plot_dot(res_pos, "Positive selection", "positive");       plot_bar(res_pos, "Positive selection", "positive") }
if (!is.null(res_kle2))      { plot_dot(res_kle2, "Ks <= 2", "Ks_leq_2");                 plot_bar(res_kle2, "Ks <= 2", "Ks_leq_2") }
if (!is.null(res_kgt2))      { plot_dot(res_kgt2, "Ks > 2", "Ks_gt_2");                   plot_bar(res_kgt2, "Ks > 2", "Ks_gt_2") }

message("All enrichments and visualisations done!")
