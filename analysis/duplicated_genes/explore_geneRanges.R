library(GenomicRanges)

# -- count nas in df$chromosome
sum(is.na(dup_df$chromosome))
sum(is.na(dup_df$start_pos))
sum(is.na(dup_df$end_pos))
# -- all are 460, removing these nas
dup_full_df <- dup_df[!is.na(dup_df$chromosome) & !is.na(dup_df$start_pos) & !is.na(dup_df$end_pos), ]

GR_duplicated <- GRanges(
  seqnames = dup_full_df$chromosome,
  ranges = IRanges(start = dup_full_df$start_pos, end = dup_full_df$end_pos),
  gene_id = dup_full_df$peptide_id,
  strand = dup_full_df$strand
)
GR_duplicated

# make gr ranges for prot_df
prot_full_df <- prot_df[!is.na(prot_df$chromosome) & !is.na(prot_df$start_pos) & !is.na(prot_df$end_pos), ]
GR_proteins <- GRanges(
  seqnames = prot_full_df$chromosome,
  ranges = IRanges(start = prot_full_df$start_pos, end = prot_full_df$end_pos),
  gene_id = prot_full_df$peptide_id,
  strand = prot_full_df$strand
)
GR_proteins

# -- make gff3
header='Chr\tGene\tStart\tEnd'

# Create a minimal MCScanX-compatible GFF
gff_mcscanx <- data.frame(
  chr = seqnames(GR_duplicated),
  start = start(GR_duplicated),
  end = end(GR_duplicated),
  gene_id = mcols(GR_duplicated)$gene_id
)

# Write to file
colnames(gff_mcscanx) <- c("Chr", "Gene", "Start", "End")
write.table(
  gff_mcscanx,
  file = "files/duplicated.gff",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# -- same for proteins
gff_proteins <- data.frame(
  chr = seqnames(GR_proteins),
  start = start(GR_proteins),
  end = end(GR_proteins),
  gene_id = mcols(GR_proteins)$gene_id
)

# Write to file
colnames(gff_proteins) <- c("Chr", "Gene", "Start", "End")
write.table(
  gff_proteins,
  file = "files/proteins.gff",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# -- circos prep


# Convert GRanges to DataFrame
df <- as.data.frame(GR_duplicated)

# ----- 1. CREATE KARYOTYPE FILE -----
karyo <- df %>%
  group_by(seqnames) %>%
  summarise(chr_start = min(start), chr_end = max(end)) %>%
  mutate(chr_id = row_number(), color = "black")

karyotype_lines <- paste(
  "chr -",
  karyo$seqnames,
  karyo$chr_id,
  karyo$chr_start,
  karyo$chr_end,
  karyo$color
)

write.table(
  karyotype_lines,
  file = "files/karyotype.txt",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

# ----- 2. EXPORT DATA VALUES FOR CIRCOS PLOT -----
# Pick a numeric column automatically if you have one
numeric_cols <- names(df)[sapply(df, is.numeric)]
numeric_cols <- setdiff(numeric_cols, c("start", "end", "width"))  # ignore structural columns

if (length(numeric_cols) == 0) {
  df$value <- 1   # fallback value if no numeric column exists
  value_col <- "value"
} else {
  value_col <- numeric_cols[1]   # use the first numeric column
}

circos_data <- df[, c("seqnames", "start", "end", value_col)]

write.table(
  circos_data,
  file = "files/data.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

# ----- 3. OUTPUT INFO -----
cat("Files created:\n")
cat(" - files/karyotype.txt\n")
cat(" - files/data.txt (for histogram/scatter/heatmap)\n")
cat("Done.\n")
