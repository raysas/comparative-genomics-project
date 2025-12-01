library(GenomicRanges)
library(dplyr)

# -- make a df of TAG pairs (geneA geneB family tag)
# -- fir
# # -- then make pairs

merged_df <- dup_full_df %>%
  inner_join(families_df, by = c("peptide_id"="geneName")) %>%
  arrange(chromosome, start_pos)

gr_all <- GRanges(
  seqnames = merged_df$chromosome,
  ranges = IRanges(start = merged_df$start_pos, end = merged_df$end_pos),
  strand = merged_df$strand
)

mcols(gr_all)$peptide_id <- merged_df$peptide_id
mcols(gr_all)$family<- merged_df$family

all_tag_results <- list()

# -------------------------------
# 3. TAG detection (global)
# -------------------------------

detect_TAGs_global <- function(gr, family, max_spacers = 5) {
  cat("running TAG detection for family:", family, "with max_spacers =", max_spacers, "\n")
  
  # Convert GRanges → sorted dataframe
  df <- data.frame(
    seqnames = as.character(seqnames(gr)),
    start = start(gr),
    end = end(gr),
    peptide_id = mcols(gr)$peptide_id,
    family = mcols(gr)$family,
    stringsAsFactors = FALSE
  )
  cat("GRanges to df")
  
  df <- df[order(df$seqnames, df$start), ]
  cat(" - sorted df\n")
  
  # Find genes of this family
  family_hits <- which(df$family == family)
  if (length(family_hits) == 0) return(NULL)
  
  tag_id <- 0
  prev_hit <- NULL
  df$TAG <- NA
  
  cat("Detecting TAGs...\n")
  
  for (i in family_hits) {
    if (is.null(prev_hit)) {
      # Start first TAG
      tag_id <- tag_id + 1
      df$TAG[i] <- tag_id
      prev_hit <- i
    } else {
      # Count number of genes between previous and current
      spacers <- i - prev_hit - 1
      if (spacers <= max_spacers) {
        # Same TAG
        df$TAG[i] <- tag_id
      } else {
        # New TAG
        tag_id <- tag_id + 1
        df$TAG[i] <- tag_id
      }
      prev_hit <- i
    }
  }
  cat('finished with family hit iterations')
  
  # Return only family genes
  cat("Returning results for family:", family, "\n")
  df[df$family == family, ]
}

# -------------------------------
# 4. Run TAG detection for all families
# -------------------------------

all_families <- unique(mcols(gr_all)$family)
tag_results_list <- lapply(all_families, function(fam) {
  detect_TAGs_global(gr_all, fam, max_spacers = 5)
})
tag_results_list <- tag_results_list[!sapply(tag_results_list, is.null)]
tag_results <- bind_rows(tag_results_list)

# -------------------------------
# 5. Flag singletons and optionally set TAG = 0
# -------------------------------

tag_sizes <- tag_results %>%
  group_by(family, TAG) %>%
  summarise(n_genes = n(), .groups = "drop")
tag_results <- tag_results %>%
  left_join(tag_sizes, by = c("family", "TAG")) %>%
  mutate(TAG = ifelse(n_genes == 1, 0, TAG)) %>%
  select(-n_genes)

# -------------------------------
# 6. Final table
# -------------------------------

tag_results <- tag_results %>%
  arrange(seqnames, start)
head(tag_results)

# -- count how many TAGs we have
num_TAGs <- sum(tag_results$TAG != 0)
cat("Number of TAGs detected ot of total:", num_TAGs,"/",nrow(tag_results), "\n")

all_tag_results[["spacer_5"]] <- tag_results

for (spacer in 0:10) {
  if (spacer == 5) next  # already done
  
  cat("Processing max_spacer =", spacer, "\n")
  
  tag_results_list <- lapply(all_families, function(fam) {
    detect_TAGs_global(gr_all, fam, max_spacers = spacer)
  })
  
  tag_results_list <- tag_results_list[!sapply(tag_results_list, is.null)]
  tag_results <- bind_rows(tag_results_list)
  
  # Flag singletons
  tag_sizes <- tag_results %>%
    group_by(family, TAG) %>%
    summarise(n_genes = n(), .groups = "drop")
  
  tag_results <- tag_results %>%
    left_join(tag_sizes, by = c("family", "TAG")) %>%
    mutate(TAG = ifelse(n_genes == 1, 0, TAG)) %>%
    select(-n_genes) %>%
    arrange(seqnames, start)
  
  all_tag_results[[paste0("spacer_", spacer)]] <- tag_results
  
  cat("for spacer=",spacer,": Number of TAGs detected ot of total:", sum(tag_results$TAG != 0),"/",nrow(tag_results), "\n")
}




#------------------------------------------------------------------------------------
#--------------------------------- Experimenting ------------------------------------
#------------------------------------------------------------------------------------


run_TAG_detection <- function(gr_all, max_spacers = 5) {
  # 1. Collect families
  cat("step 1: Collecting families...\n")
  all_families <- unique(mcols(gr_all)$family)

  # 2. Run detection for each family
  cat("step 2: Running TAG detection for max_spacers =", max_spacers, "...\n")
  tag_results_list <- lapply(all_families, function(fam) {
    detect_TAGs_global(gr_all, fam, max_spacers = max_spacers)
  })
  
  # 3. Remove NULL results and bind
  cat("step 3: Combining results...\n")
  tag_results_list <- tag_results_list[!sapply(tag_results_list, is.null)]
  tag_results <- bind_rows(tag_results_list)
  cat("Total genes with TAGs detected:", nrow(tag_results), "\n")
  
  # 4. Flag singletons (TAG = 0 if only one gene in family)
  cat("step 4: Flagging singletons...\n")
  tag_sizes <- tag_results %>%
    group_by(family, TAG) %>%
    summarise(n_genes = n(), .groups = "drop")
  
  tag_results <- tag_results %>%
    left_join(tag_sizes, by = c("family", "TAG")) %>%
    mutate(TAG = ifelse(n_genes == 1, 0, TAG)) %>%
    select(-n_genes)
  
  # 5. Final table sorted
  cat("step 5: Sorting final results...\n")
  tag_results <- tag_results %>%
    arrange(seqnames, start)
  
  # 6. Report summary
  cat("step 6: Summary...\n")
  num_TAGs <- sum(tag_results$TAG != 0)
  cat("Number of TAGs detected out of total:", num_TAGs, "/", nrow(tag_results), "\n")
  
  return(tag_results)
}
res0 <- run_TAG_detection(gr_all, max_spacers = 0)
res1 <- run_TAG_detection(gr_all, max_spacers = 1)
res2 <- run_TAG_detection(gr_all, max_spacers = 2)
res3 <- run_TAG_detection(gr_all, max_spacers = 3)
res4 <- run_TAG_detection(gr_all, max_spacers = 4)
res5 <- run_TAG_detection(gr_all, max_spacers = 5)
res6 <- run_TAG_detection(gr_all, max_spacers = 6)
res7 <- run_TAG_detection(gr_all, max_spacers = 7)
res8 <- run_TAG_detection(gr_all, max_spacers = 8)
res9 <- run_TAG_detection(gr_all, max_spacers = 9)
res10 <- run_TAG_detection(gr_all, max_spacers = 10)

save(res0, res1, res2, res3, res4, res5, res6, res7, res8, res9, res10,
     file = "files/tag_spacers_experiment.Rdata")















# -- looping over all spacers
library(GenomicRanges)
library(dplyr)

# Prepare merged GRanges object as you already have
merged_df <- dup_full_df %>%
  inner_join(families_df, by = c("peptide_id"="geneName")) %>%
  arrange(chromosome, start_pos)

gr_all <- GRanges(
  seqnames = merged_df$chromosome,
  ranges = IRanges(start = merged_df$start_pos, end = merged_df$end_pos),
  strand = merged_df$strand
)

mcols(gr_all)$peptide_id <- merged_df$peptide_id
mcols(gr_all)$family <- merged_df$family

# -------------------------------
# Run TAG detection for spacers 0 to 10
# -------------------------------

max_spacer_values <- 0:10
all_tag_results <- list()

for (spacer in max_spacer_values) {
  cat("Processing max_spacer =", spacer, "\n")
  
  tag_results_list <- lapply(all_families, function(fam) {
    detect_TAGs_global(gr_all, fam, max_spacers = spacer)
  })
  
  # Remove empty results
  tag_results_list <- tag_results_list[!sapply(tag_results_list, is.null)]
  
  # Combine all families
  tag_results <- bind_rows(tag_results_list)
  
  # Flag singletons
  tag_sizes <- tag_results %>%
    group_by(family, TAG) %>%
    summarise(n_genes = n(), .groups = "drop")
  
  tag_results <- tag_results %>%
    left_join(tag_sizes, by = c("family", "TAG")) %>%
    mutate(TAG = ifelse(n_genes == 1, 0, TAG)) %>%
    select(-n_genes) %>%
    arrange(seqnames, start)
  
  # Store in list
  all_tag_results[[paste0("spacer_", spacer)]] <- tag_results
  
  cat("for spacer=",spacer,": Number of TAGs detected ot of total:", sum(tag_results$TAG != 0),"/",nrow(tag_results), "\n")
}

# -------------------------------
# Save all results in one RData file
# -------------------------------

save(all_tag_results, file = "files/tag_spacers_experiment.Rdata")
###################################################################################3
#######################################################################################
######################################################################################
# -- parallelized versio --
library(dplyr)
library(GenomicRanges)
library(parallel) # for mclapply
library(dplyr)

library(future.apply)

plan(multisession)   # or multicore on Linux/macOS

all_tag_results <- vector("list", length(max_spacer_values))
names(all_tag_results) <- paste0("spacer_", max_spacer_values)

max_spacer_values <- 0:10
for (spacer in max_spacer_values) {
  cat("Processing max_spacer =", spacer, "\n")
  
  # PARALLEL apply
  tag_results_list <- future_lapply(all_families, function(fam) {
    detect_TAGs_global(gr_all, fam, max_spacers = spacer)
  })
  
  tag_results_list <- tag_results_list[!sapply(tag_results_list, is.null)]
  tag_results <- bind_rows(tag_results_list)
  
  # Compute TAG sizes once
  tag_results <- tag_results %>%
    group_by(family, TAG) %>%
    mutate(TAG = ifelse(n() == 1, 0, TAG)) %>%
    ungroup() %>%
    arrange(seqnames, start)
  
  all_tag_results[[paste0("spacer_", spacer)]] <- tag_results
  
  cat(
    "for spacer =", spacer, 
    ": Number of TAGs detected out of total:", 
    sum(tag_results$TAG != 0), "/", nrow(tag_results), "\n"
  )
}
# Save all results
save(all_tag_results, file = "files/tag_spacers_experiment.Rdata")


##############################################


