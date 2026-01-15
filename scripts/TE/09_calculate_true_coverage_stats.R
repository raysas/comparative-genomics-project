# ------------------------------------------------------------------
# Script Name: 09_calculate_true_coverage_stats.R
#
# Description: 
#   PART 1: DATA EXTRACTION
#   Calculates TRUE genomic coverage metrics for Transposable Elements.
#   - METHOD: Uses `GenomicRanges` to calculate both:
#       1. Simple Sum (Total length of all annotations, including overlaps).
#       2. Reduced Sum (Physical DNA footprint, overlaps merged).
#   - METRICS: 
#       - Genome_Pct: % of physical genome occupied.
#       - Nesting_Index: Ratio of Simple Sum / True BP (Measure of stacking).
#   - INPUT: GFF3 (TEs) + TSV (Chromosome Lengths).
#   - OUTPUT: A comprehensive CSV table ("TE_true_coverage_stats.csv").
# ------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(scales)

# ==================================================================
# 1. LOAD DATA & GENOME SIZE
# ==================================================================
gff_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/TEAnnotationFinal.gff3"
len_path <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/data/chr_lengths.tsv"

# A. Calculate Total Genome Size (Denominator)
# We use the sum of the 20 chromosomes (the analyzed space) as the "True Genome Size"
chr_ref <- read.table(len_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
TOTAL_GENOME_BP <- sum(as.numeric(chr_ref[,2])) 

print(paste("Total Genome Size (Chr 1-20):", comma(TOTAL_GENOME_BP), "bp"))

# B. Load TE Data
gff_data <- read.table(gff_path, sep="\t", comment.char="#", quote="", stringsAsFactors=FALSE)
raw_data <- gff_data %>% dplyr::select(V1, V3, V4, V5)
colnames(raw_data) <- c("Chr", "Full_Name", "Start", "End")

# Clean and Parse
clean_data <- raw_data %>%
  tidyr::separate(Full_Name, into = c("Class_Raw", "Order", "Superfamily"), sep = "/", fill = "right") %>%
  dplyr::mutate(
    Class_Group = case_when(
      grepl("Class II", Class_Raw) ~ "Class II (DNA)",
      grepl("Class I", Class_Raw) ~ "Class I (Retro)",
      TRUE ~ "Unknown"
    ),
    Order = ifelse(Order == "-" | is.na(Order) | Order == "", "Unknown", Order),
    # Clean Superfamily name
    Superfamily = ifelse(Superfamily == "-" | is.na(Superfamily) | Superfamily == "", "Unspecified", Superfamily),
    # Create Full Label
    Type_Label = paste0(Order, "/", Superfamily)
  ) %>%
  dplyr::filter(!grepl("scaffold|contig|chl|mito|Unknown", Chr, ignore.case = TRUE)) %>%
  dplyr::filter(Start > 0 & End > Start)

# ==================================================================
# 2. CALCULATE TRUE COVERAGE (Using GenomicRanges)
# ==================================================================

# Helper function to calculate BOTH simple and reduced widths
calc_coverage_stats <- function(chr, start, end) {
  if (length(chr) == 0) return(data.frame(Simple_Sum_BP=0, Reduced_BP=0))
  
  # Create GRanges
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  
  # Simple Sum (Total annotation length, including overlaps)
  simple_sum <- sum(as.numeric(width(gr)))
  
  # Reduced Sum (Physical footprint, overlaps merged)
  reduced_sum <- sum(as.numeric(width(reduce(gr))))
  
  return(data.frame(Simple_Sum_BP = simple_sum, Reduced_BP = reduced_sum))
}

print("Calculating Global Stats...")
# --- A. Global TE Coverage ---
global_stats <- calc_coverage_stats(clean_data$Chr, clean_data$Start, clean_data$End)
total_te_bp <- global_stats$Reduced_BP

print("--- GLOBAL STATISTICS ---")
print(paste("Simple Sum (Nested/Redundant):", comma(global_stats$Simple_Sum_BP), "bp"))
print(paste("Reduced Sum (Physical Footprint):", comma(global_stats$Reduced_BP), "bp"))
print(paste("Physical Genome Coverage:", round((global_stats$Reduced_BP / TOTAL_GENOME_BP) * 100, 2), "%"))
print(paste("Redundant Content Ratio:", round((global_stats$Simple_Sum_BP / TOTAL_GENOME_BP) * 100, 2), "%"))

print("Calculating Class Stats...")
# --- B. Per Class Coverage ---
class_stats <- clean_data %>%
  dplyr::group_by(Class_Group) %>%
  dplyr::do(calc_coverage_stats(.$Chr, .$Start, .$End)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Category = "Class") %>%
  dplyr::rename(Type = Class_Group, True_BP = Reduced_BP)

print("Calculating Superfamily Stats...")
# --- C. Per Superfamily Coverage ---
# Group by Class_Group, Order AND Superfamily
sf_stats <- clean_data %>%
  dplyr::group_by(Class_Group, Order, Superfamily, Type_Label) %>%
  dplyr::do(calc_coverage_stats(.$Chr, .$Start, .$End)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Category = "Superfamily") %>%
  dplyr::rename(Type = Type_Label, True_BP = Reduced_BP)

# Combine for CSV Output
total_row <- data.frame(
  Class_Group = "All", 
  Type = "Total Genome TE", 
  Simple_Sum_BP = global_stats$Simple_Sum_BP,
  True_BP = total_te_bp, 
  Category = "Total"
)

# Fix Class_Group column in class_stats before binding
class_stats$Class_Group <- class_stats$Type

# Combine
final_stats <- dplyr::bind_rows(total_row, class_stats, sf_stats)

# ==================================================================
# 3. DERIVE METRICS & SAVE
# ==================================================================
# Calculate Percentage of Genome
final_stats$Genome_Pct <- round((final_stats$True_BP / TOTAL_GENOME_BP) * 100, 3) 
final_stats$Mb <- round(final_stats$True_BP / 1000000, 3)

# Calculate Nesting Index
# Formula: Simple Sum / True BP
# 1.0 = No Nesting. 1.5 = 50% "extra" nested sequence.
final_stats$Nesting_Index <- round(final_stats$Simple_Sum_BP / final_stats$True_BP, 2)

# Save Table
out_csv <- "C:/Users/dell/OneDrive/Master/M2/Comparative_Genomics/comparative-genomics-project/output/TE_true_coverage_stats.csv"
write.csv(final_stats, out_csv, row.names=FALSE)
print(paste("Data Extraction Complete. Stats saved to:", out_csv))