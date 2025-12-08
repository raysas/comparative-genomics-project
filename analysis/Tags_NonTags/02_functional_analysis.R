# Input: PANTHER GO-Slim enrichment tables
# Output: GOslim_low_dotplot.png, GOslim_high_dotplot.png

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(patchwork)
library(cowplot)

# 1. Read PANTHER GO file
read_panther_file <- function(path) {

  lines <- readLines(path, encoding = "UTF-8")
  header_line <- grep("PANTHER GO-Slim", lines)[1]
  if (is.na(header_line)) stop(paste("No header found in", path))

  df <- read.delim(path, header = TRUE, sep = "\t", skip = header_line - 1)

  colnames(df) <- c("Term", "RefList", "QueryList", "Expected",
                    "OverUnder", "FoldEnrichment", "Pvalue", "FDR")

  df <- df %>% filter(str_detect(Term, "GO:"))

  df$GO_ID   <- str_extract(df$Term, "GO:\\d+")
  df$GO_Name <- str_replace(df$Term, " \\(GO:\\d+\\)", "")
  df$GeneRatio <- df$QueryList / df$RefList

  df %>% drop_na(FDR, GeneRatio)
}

# 2. Select Top 10 enriched terms (by FDR)
top10 <- function(df) {
  df %>%
    arrange(FDR) %>%
    slice_head(n = 10) %>%
    arrange(GeneRatio)
}

# 3. Dotplot (enrichplot-style)
dotplot_panther <- function(df, title) {

  df$GO_Name <- factor(df$GO_Name, levels = df$GO_Name)
  max_count <- max(df$QueryList)
  size_range <- if (max_count > 3000) c(1.5, 6) else c(2.5, 8)

  p <- ggplot(df,
              aes(x = GeneRatio,
                  y = GO_Name,
                  size = QueryList,
                  color = -log10(FDR))) +

    geom_point(alpha = 0.9) +

    scale_color_gradientn(
      colours = c("#9ECAE1", "#084594", "#7A0177"),
      name = "-log10(p.adjust)"
    ) +

    scale_size_continuous(
      name = "Count",
      breaks = pretty(df$QueryList, 3),
      range = size_range
    ) +

    guides(
      color = guide_colorbar(order = 1),
      size  = guide_legend(order = 2)
    ) +

    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title.y = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      legend.spacing.y = unit(0.3, "cm"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9)
    ) +

    labs(x = "GeneRatio", title = title)

  return(p)
}

# 4. Build 9-panel figure for a given stringency
plot_stringency <- function(stringency) {

  categories <- c("BP", "CC", "MF")
  families   <- c("singletons", "TAGs", "non_TAGs")  

  input_dir  <- "../../output/Tags_NonTags/GOenrich/"
  output_dir <- "../../figures/Tags_NonTags/"
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  plots <- list()
  k <- 1

  for (cat in categories) {
    for (fam in families) {

      file_path <- file.path(
        input_dir,
        paste0(cat, "_", fam, "_", stringency, ".txt")
      )

      df  <- read_panther_file(file_path)
      df2 <- top10(df)

      cat_label <- case_when(
        cat == "BP" ~ "Biological process",
        cat == "CC" ~ "Cellular component",
        cat == "MF" ~ "Molecular function"
      )

      fam_label <- case_when(
        fam == "TAGs"    ~ "TAGs",
        fam == "singletons" ~ "Singletons",
        fam == "non_TAGs"    ~ "Non_TAGs"
      )

      title <- paste(cat_label, "–", fam_label)

      p <- dotplot_panther(df2, title)
      p <- p + theme(plot.margin = margin(10, 10, 10, 10))

      plots[[k]] <- p
      k <- k + 1
    }
  }

  aligned <- align_plots(
      plots[[1]], plots[[2]], plots[[3]],
      plots[[4]], plots[[5]], plots[[6]],
      plots[[7]], plots[[8]], plots[[9]],
      align = "vh",
      axis = "tblr"
  )

  final_plot <- plot_grid(
      aligned[[1]], aligned[[2]], aligned[[3]],
      aligned[[4]], aligned[[5]], aligned[[6]],
      aligned[[7]], aligned[[8]], aligned[[9]],
      ncol = 3,
      rel_heights = c(1, 1, 1),
      rel_widths = c(1, 1, 1)
  )

  final_plot <- final_plot +
    plot_annotation(
      title = paste(str_to_title(stringency), "stringency"),
      theme = theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5))
    )

  outfile <- file.path(output_dir,
                       paste0("GOslim_", stringency, "_dotplot.png"))

  ggsave(outfile, final_plot, width = 30, height = 14, dpi = 300)
  message("Saved: ", outfile)
}

# 5. Run
plot_stringency("low")
plot_stringency("high")
