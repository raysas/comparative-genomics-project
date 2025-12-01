#!/usr/bin/env Rscript

# -- script: proximal_genes_by_family.R
# -- purpose: for each gene family, find all gene pairs proximal by n kilobases

suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(readr)
  library(argparse)
  library(parallel)
})

# -- argument parser
parser <- ArgumentParser(description = "find proximal gene pairs by family")
parser$add_argument("-g", "--genes", default = "output/statistics/duplicated_genes_info.csv",
                    help = "gene annotation csv file (default: %(default)s)")
parser$add_argument("-f", "--family", default = "output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv",
                    help = "family mapping csv file with geneName and family_id (default: %(default)s)")
parser$add_argument("-d", "--distance", type = "double", default = 500,
                    help = "maximum distance in kilobases for proximity (default: %(default)s kb)")
parser$add_argument("-o", "--output", default = "output/duplication_classes/proximal/proximal_%dkb.csv",
                                        help = "output csv file, %d will be replaced by distance (default: %(default)s)")

args <- parser$parse_args()

genes_file <- args$genes
family_file <- args$family
proximal_kb <- args$distance
output_file <- sprintf(args$output, proximal_kb)

# -- load data
genes_df <- read.csv(genes_file)
family_df <- read.table(family_file, header = TRUE)

print(paste("-- total genes loaded:", nrow(genes_df)))

# -- merge genes with family info
genes_df <- genes_df %>%
  inner_join(family_df, by = c("peptide_id" = "geneName"))

# -- remove NA chromosome from genes
genes_df <- genes_df %>%
  filter(!is.na(chromosome))
print(paste("-- total genes after removing NA chromosome:", nrow(genes_df)))

# -- create granges object
gr <- GRanges(
  seqnames = genes_df$chromosome,
  ranges = IRanges(start = genes_df$start_pos, end = genes_df$end_pos),
  strand = ifelse(genes_df$strand == 1, "+", "-"),
  gene_id = genes_df$gene_id,
  peptide_id = genes_df$peptide_id,
  family_id = genes_df$family
)

print(paste("-- created GRanges with total genes:", length(gr)))

print("-- combining proximal pairs from all families...")

# -- function to get proximal pairs per family
get_proximal_pairs <- function(family_df, gr, max_dist) {

    ids <- family_df$peptide_id
    if (length(ids) < 2) return(list())   # -- must return list

    # -- get the relevant GRanges entries
    gr_sub <- gr[gr$peptide_id %in% ids]

    if (length(gr_sub) < 2) return(list())

    combs <- combn(gr_sub$peptide_id, 2, simplify = FALSE)

    res <- lapply(combs, function(x) {
        g1 <- gr_sub[gr_sub$peptide_id == x[1]]
        g2 <- gr_sub[gr_sub$peptide_id == x[2]]

        # -- only consider same chromosome
        if (as.character(seqnames(g1)) != as.character(seqnames(g2)))
            return(NULL)

        dist <- abs(start(g1) - start(g2))

        if (dist <= max_dist) {
            return(data.frame(
                geneA = x[1],
                geneB = x[2],
                chromosome = as.character(seqnames(g1)),
                distance = dist
            ))
        } else {
            return(NULL)
        }
    })

    # -- remove NULL entries
    res <- res[!sapply(res, is.null)]
    if (length(res) == 0) return(list())
    return(res)
}


# -- loop over families
# proximal_pairs_list <- lapply(unique(gr$family_id), function(fam) {
#   gr_fam <- gr[gr$family_id == fam]
#   get_proximal_pairs(gr_fam, gr, max_dist = proximal_kb * 1000)
#   print(paste("-- processed family:", fam))
# })

# -- parallel version

num_cores <- 8

proximal_pairs_list <- mclapply(
  unique(gr$family_id),
  function(fam) {
    gr_fam <- gr[gr$family_id == fam]
    res <- get_proximal_pairs(gr_fam, gr, max_dist = proximal_kb * 1000)
    message(paste("-- processed family:", fam))
    return(res)
  },
  mc.cores = num_cores
)


print("-- done! making final dataframe...")

proximal_pairs_df <- do.call(rbind, proximal_pairs_list)

# -- save output
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.csv(proximal_pairs_df, output_file, row.names = FALSE)
cat("saved proximal gene pairs to:", output_file, "\n")
