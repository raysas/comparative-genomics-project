library(dplyr)
library(readr)
library(gprofiler2)
library(ggplot2)
library(forcats)

# ==========================================================
# 1. LOAD DATA (fixed paths)
# ==========================================================

df <- read_delim(
  "../../output/ks_results_high/ks_with_GLYMA_ID_GO_high.tsv",
  delim = "\t",
  col_types = cols()
)

out_enrich_dir <- "../../output/enrichment_gprofiler_high"
out_fig_dir    <- "../../figures/enrichment_gprofiler_high"

dir.create(out_enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir,    recursive = TRUE, showWarnings = FALSE)

df <- df %>% 
  filter(!is.na(gene)) %>% 
  mutate(GLYMA_ID = as.character(GLYMA_ID))


# ==========================================================
# 2. GENE SETS
# ==========================================================

recent_genes     <- df %>% filter(age_class == "Recent") %>% pull(GLYMA_ID) %>% unique()
old_genes        <- df %>% filter(age_class == "Old")    %>% pull(GLYMA_ID) %>% unique()
purifying_genes  <- df %>% filter(selection_class == "Purifying_selection") %>% pull(GLYMA_ID) %>% unique()
positive_genes   <- df %>% filter(selection_class == "Positive_selection") %>% pull(GLYMA_ID) %>% unique()

df <- df %>% mutate(ks_group = ifelse(ks_used <= 2, "<=2", ">2"))
ks_le2_genes <- df %>% filter(ks_group == "<=2") %>% pull(GLYMA_ID) %>% unique()
ks_gt2_genes <- df %>% filter(ks_group == ">2")  %>% pull(GLYMA_ID) %>% unique()


# ==========================================================
# 3. RUN g:PROFILER
# ==========================================================

run_gprof <- function(genes, name){
  
  if (length(genes) == 0){
    message("Skipping ", name, ": empty list")
    return(NULL)
  }
  
  gp <- gost(
    query    = genes,
    organism = "gmax",
    sources  = c("GO:BP", "GO:MF", "GO:CC")
  )
  
  if (is.null(gp$result)){
    message("No enrichment for ", name)
    return(NULL)
  }

  sanitize_df <- function(x){
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    for (colname in names(df)){
      col_vals <- df[[colname]]
      is_problem <- is.list(col_vals) || any(sapply(col_vals,
                        function(z) is.list(z) || length(z) > 1))
      if (is_problem){
        df[[colname]] <- sapply(col_vals, function(z){
          if (is.null(z)) return(NA_character_)
          paste(unlist(z), collapse = ";")
        })
      } else {
        df[[colname]] <- as.character(col_vals)
      }
    }
    df
  }
  
  res_df <- sanitize_df(gp$result)
  write.csv(res_df, file.path(out_enrich_dir, paste0(name, "_gprof.csv")), row.names = FALSE)
  
  return(res_df)
}

res_recent <- run_gprof(recent_genes,    "Recent")
res_old    <- run_gprof(old_genes,       "Old")
res_pur    <- run_gprof(purifying_genes, "Purifying")
res_pos    <- run_gprof(positive_genes,  "Positive")
res_kle2   <- run_gprof(ks_le2_genes,    "Ks_leq_2")
res_kgt2   <- run_gprof(ks_gt2_genes,    "Ks_gt_2")


# ==========================================================
# 4. VISUALIZATION FUNCTIONS
# ==========================================================

# ---- A. Dotplot (BP only) ----
plot_dot <- function(res, title, filename){
  df <- res %>% 
    filter(source == "GO:BP") %>% 
    mutate(p_value = as.numeric(p_value)) %>% 
    filter(!is.na(p_value)) %>% 
    arrange(p_value) %>% 
    slice_head(n = 10) %>% 
    mutate(term_name = fct_reorder(term_name, -log10(p_value)))
  
  p <- ggplot(df, aes(x = term_name, y = -log10(p_value))) +
    geom_point(size = 4, color = "#084594") +
    coord_flip() +
    theme_minimal() +
    labs(title = title, x = "GO term (BP)", y = "-log10(p-value)")
  
  ggsave(file.path(out_fig_dir, paste0(filename, "_dotplot.png")), p, width = 10, height = 6)
}

# ---- B. Barplot (BP only) ----
plot_bar <- function(res, title, filename){
  df <- res %>% 
    filter(source == "GO:BP") %>% 
    mutate(p_value = as.numeric(p_value)) %>% 
    filter(!is.na(p_value)) %>% 
    arrange(p_value) %>% 
    slice_head(n = 10) %>% 
    mutate(term_name = fct_reorder(term_name, -log10(p_value)))
  
  p <- ggplot(df, aes(x = term_name, y = -log10(p_value))) +
    geom_bar(stat = "identity", fill = "#7A0177") +
    coord_flip() +
    theme_minimal() +
    labs(title = title, x = "GO term (BP)", y = "-log10(p-value)")
  
  ggsave(file.path(out_fig_dir, paste0(filename, "_barplot.png")), p, width = 10, height = 6)
}

# ---- C. Bubble Plot (BP + MF) ----
plot_bubble <- function(res, title, filename, ont = "GO:BP"){
  df <- res %>% 
    filter(source == ont) %>% 
    mutate(
      p_value           = as.numeric(p_value),
      intersection_size = as.numeric(intersection_size),
      query_size        = as.numeric(query_size)
    ) %>% 
    filter(!is.na(p_value), p_value > 0,
           !is.na(intersection_size), !is.na(query_size)) %>%
    mutate(
      GeneRatio = intersection_size / query_size,
      neglog10p = -log10(p_value)
    ) %>% 
    arrange(p_value) %>% 
    slice_head(n = 10) %>% 
    mutate(term_name = fct_reorder(term_name, GeneRatio))
  
  p <- ggplot(df, aes(
    x = GeneRatio,
    y = term_name,
    size = intersection_size,
    colour = neglog10p
  )) +
    geom_point(alpha = 0.85) +
    scale_colour_gradient(low = "#81199b", high = "#401191") +
    theme_minimal() +
    labs(
      title = title,
      x = "GeneRatio",
      y = "",
      colour = "-log10(p-value)",
      size = "Gene count"
    )
  
  suffix <- ifelse(ont == "GO:BP", "BP", "MF")

  ggsave(file.path(out_fig_dir, paste0(filename, "_", suffix, "_bubble.png")),
      p, width = 10, height = 6)

}


# ==========================================================
# 5. GENERATE ALL PLOTS (BP + MF bubble only)
# ==========================================================

make_plots <- function(res, title, filename){
  if (!is.null(res)){

    # BP plots
    plot_dot(res, title, filename)
    plot_bar(res, title, filename)
    plot_bubble(res, title, filename, ont = "GO:BP")

    # MF bubble only
    plot_bubble(res, paste(title, "(MF)"), filename, ont = "GO:MF")

  }
}

make_plots(res_recent, "Recent duplicates",       "recent")
make_plots(res_old,    "Old duplicates",          "old")
make_plots(res_pur,    "Purifying selection",     "purifying")
make_plots(res_pos,    "Positive selection",      "positive")


message("All enrichments and visualisations done!")
