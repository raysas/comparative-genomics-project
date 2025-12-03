#!/usr/bin/env Rscript

# ---------------------------
# -- script: tag_pair_orientation.R
# -- purpose: Generate all tag-based gene pairs, compute their orientation,
#           and summarize orientation counts.
# Rscript scripts/TAGs_pairs_orientation.R --stringency low
# Rscript scripts/TAGs_pairs_orientation.R --stringency high
# ---------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(argparse)
})

# ---------------------------
# -- argument parser
# ---------------------------
parser <- ArgumentParser(description = "Compute TAG gene pairs with orientation")
parser$add_argument("-t", "--tags", default = "output/duplication_classes/TAGs/TAGs_1.tsv",
                    help = "Input TAGs file (default: %(default)s)")
parser$add_argument("-d", "--dup", default = "output/info/duplicated_genes_info.csv",
                    help = "Input duplication full data file (default: %(default)s)")
parser$add_argument("-o", "--output", default = "output/duplication_classes/TAGs/TAG_gene_pairs.tsv",
                    help = "Output file for TAG gene pairs (default: %(default)s)")
parser$add_argument("--stringency", choices = c("low", "high"), default = NULL,
                    help = "-- stringency level for duplicated genes and families [default: %(default)s]")

args <- parser$parse_args()

tags_input_file <- args$tags
dup_input_file <- args$dup
output_file <- args$output

if (!is.null(args$stringency)) {
  if (args$stringency == "low") {
    print("-- using low stringency duplicated genes, overriding any provided file paths: id30, cov50, evalue1e-10")
    dup_input_file <- "output/info/duplicated_genes_info_id30_qcov50_scov50_evalue1e-10_wcol12.csv"
    tags_input_file <- "output/duplication_classes/TAGs/low/TAGs_1.tsv"
    output_file <- "output/duplication_classes/TAGs/low/TAG_gene_pairs.tsv"
  } else if (args$stringency == "high") {
    print("-- using high stringency duplicated genes, overriding any provided file paths: filter id50, cov70, evalue1e-10")
      dup_input_file <- "output/info/duplicated_genes_info_id50_qcov70_scov70_evalue1e-10_wcol12.csv"
      tags_input_file <- "output/duplication_classes/TAGs/high/TAGs_1.tsv"
      output_file <- "output/duplication_classes/TAGs/high/TAG_gene_pairs.tsv"
  }
}

# ---------------------------
# -- load data
# ---------------------------
tags_df <- read.table(tags_input_file, header = TRUE, stringsAsFactors = FALSE)
dup_full_df <- read.csv(dup_input_file, header = TRUE, stringsAsFactors = FALSE)
# dup_full_df <- dup_full_df %>%
#   filter(!is.na(start_pos) & !is.na(end_pos) & !is.na(strand))

# print(colnames(dup_full_df))

# ---------------------------
# Create tag IDs
# ---------------------------
tags_df$tag_id <- ifelse(tags_df$TAG == 0,
                         0,
                         paste0(tags_df$family, "_TAG", tags_df$TAG))

print(paste0("-- total TAG genes identified: ", sum(tags_df$TAG != 0)))
print(paste("-- total TAG arrays identified (excluding 0 in tag id):", length(unique(tags_df$tag_id[tags_df$tag_id != 0]))))

# ---------------------------
# -- generate all gene pairs per tag
# ---------------------------
tag_pairs_list <- lapply(unique(tags_df$tag_id[tags_df$tag_id != 0]), function(tid) {
  members <- tags_df$peptide_id[tags_df$tag_id == tid]
  if (length(members) < 2) return(NULL)
  combn(members, 2, FUN = function(x) {
    data.frame(
      geneA = x[1],
      geneB = x[2],
      tag_id = tid,
      startA = dup_full_df$start_pos[dup_full_df$peptide_id == x[1]],
      endA = dup_full_df$end_pos[dup_full_df$peptide_id == x[1]],
      strandA = dup_full_df$strand[dup_full_df$peptide_id == x[1]],
      startB = dup_full_df$start_pos[dup_full_df$peptide_id == x[2]],
      endB = dup_full_df$end_pos[dup_full_df$peptide_id == x[2]],
      strandB = dup_full_df$strand[dup_full_df$peptide_id == x[2]],
      stringsAsFactors = FALSE
    )
  }, simplify = FALSE)
})

tag_pairs_df <- do.call(rbind, unlist(tag_pairs_list, recursive = FALSE))

print(paste("-- total TAG gene pairs generated:", nrow(tag_pairs_df)))

# ---------------------------
# -- order pairs so that startA < startB
# ---------------------------
tag_pairs_df <- tag_pairs_df %>%
  rowwise() %>%
  mutate(
    geneA_final = ifelse(startA < startB, geneA, geneB),
    geneB_final = ifelse(startA < startB, geneB, geneA),
    startA_final = ifelse(startA < startB, startA, startB),
    endA_final = ifelse(startA < startB, endA, endB),
    strandA_final = ifelse(startA < startB, strandA, strandB),
    startB_final = ifelse(startA < startB, startB, startA),
    endB_final = ifelse(startA < startB, endB, endA),
    strandB_final = ifelse(startA < startB, strandB, strandA)
  ) %>%
  select(geneA = geneA_final,
         geneB = geneB_final,
         tag_id,
         startA = startA_final,
         endA = endA_final,
         strandA = strandA_final,
         startB = startB_final,
         endB = endB_final,
         strandB = strandB_final) %>%
  ungroup()

print(paste("-- total TAG gene pairs generated:", nrow(tag_pairs_df)))

# ---------------------------
# -- add orientation
# ---------------------------
tag_pairs_df <- tag_pairs_df %>%
  rowwise() %>%
  mutate(orientation = case_when(
    strandA == strandB ~ "tandem",
    strandA == 1 & strandB == -1 ~ "convergent",
    strandA == -1 & strandB == 1 ~ "divergent",
    TRUE ~ NA_character_
  )) %>%
  ungroup()

print(paste("-- total TAG gene pairs after ordering:", nrow(tag_pairs_df)))

# ---------------------------
# -- write output
# ---------------------------
write.table(tag_pairs_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
print(paste("-- TAG gene pairs with orientation written to:", output_file))

# ---------------------------
# -- orientation summary
# ---------------------------
orientation_counts <- tag_pairs_df %>%
  group_by(orientation) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count)) %>%
  mutate(percentage = (count / sum(count)) * 100)

print(orientation_counts)

# ---------------------------
# -- goodness of if 
# ---------------------------

observed <- table(tag_pairs_df$orientation)
expected_prop <- c(tandem = 0.5, convergent = 0.25, divergent = 0.25)
expected <- sum(observed) * expected_prop[names(observed)]

print(chisq.test(x = observed, p = expected_prop[names(observed)]))