

# -- takes arguments:
#      - stringency level: low or high (default: low)

# if (!requireNamespace(c("argparse", "dplyr"), quietly = TRUE)) {
#   install.packages(c("argparse", "dplyr"))
# }
suppressPackageStartupMessages({
    library(dplyr)
    library(argparse)
})

# -- get arguments
parser <- ArgumentParser(description = "extract singleton genes from the protein file based on duplicated genes info")
parser$add_argument("--dup_file", default = "output/info/duplicated_genes_info.csv",
                    help = "-- path to duplication dataframe CSV [default: %(default)s]")
parser$add_argument("--prot_file", default = "output/info/protein_info_longest.csv",
                    help = "-- path to proteins dataframe CSV [default: %(default)s]")
parser$add_argument("--outfile", default = "output/duplication_classes/singletons/singletons.csv",
                    help = "-- path to save the output CSV/TSV [default: output/duplication_classes/singletons/singletons.csv]")
parser$add_argument("--stringency", choices = c("low", "high"), default = NULL,
                    help = "-- stringency level for duplicated genes [default: %(default)s]")

args <- parser$parse_args()

if (is.null(args$stringency)) {
  dup_file <- args$dup_file
    prot_file <- args$prot_file
    outfile <- args$outfile
    suffix <- "default"
} else if (args$stringency == "low") {
  print("-- using low stringency duplicated genes, overriding any provided file paths: id30, cov50, evalue1e-10")
  dup_file <- "output/info/duplicated_genes_info_id30_qcov50_scov50_evalue1e-10_wcol12.csv"
  prot_file <- args$prot_file
  outfile <- paste0("output/duplication_classes/singletons/singletons_low.csv")
  suffix <- "low"
} else if (args$stringency == "high") {
  print("-- using high stringency duplicated genes, overriding any provided file paths: filter id50, cov70, evalue1e-10")
  dup_file <- "output/info/duplicated_genes_info_id50_qcov70_scov70_evalue1e-10_wcol12.csv"
  prot_file <- args$prot_file
  outfile <- paste0("output/duplication_classes/singletons/singletons_high.csv")
  suffix <- "high"
}

## ensure input files exist
if (!file.exists(dup_file)) {
  stop(paste0("Duplication input file not found: ", dup_file))
}
if (!file.exists(prot_file)) {
  stop(paste0("Proteins input file not found: ", prot_file))
}

dup_df <- tryCatch(read.csv(dup_file, stringsAsFactors = FALSE),
                   error = function(e) stop(paste0("Failed to read duplication file: ", e$message)))
dup_pep_ids<- dup_df$peptide_id

prot_df <- tryCatch(read.csv(prot_file, stringsAsFactors = FALSE),
                    error = function(e) stop(paste0("Failed to read proteins file: ", e$message)))
all_pep_ids<- prot_df$peptide_id

singleton_pep_ids<- setdiff(all_pep_ids, dup_pep_ids)
singleton_df<- prot_df[prot_df$peptide_id %in% singleton_pep_ids, ]

# create output directories if they don't exist
outdir <- dirname(outfile)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(singleton_df, file=outfile, row.names=FALSE, quote=FALSE)
print(paste0("-- saved singletons dataframe to ", outfile))

# -- now saving singletons in gene_lists/singletons/
singleton_ids_file <- paste0("output/gene_lists/singletons/singletons_", suffix, ".txt")
gene_list_dir <- dirname(singleton_ids_file)
if (!dir.exists(gene_list_dir)) {
  dir.create(gene_list_dir, recursive = TRUE, showWarnings = FALSE)
}
write.table(singleton_pep_ids, file=singleton_ids_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
print(paste0("-- saved singletons peptide IDs to ", singleton_ids_file))