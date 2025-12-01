
# ------------------ for tags df -------------------------
tags_df<-read.table('../../output/duplication_classes/TAGS/TAGs_1.tsv', header=TRUE)

tags_df$tag_id <- ifelse(tags_df$TAG == 0,
                         0,
                         paste0(tags_df$family, "_TAG", tags_df$TAG))

# --for every pair having same tag_id, make all combinations
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
      strandA= dup_full_df$strand[dup_full_df$peptide_id == x[1]],
      startB = dup_full_df$start_pos[dup_full_df$peptide_id == x[2]],
      endB = dup_full_df$end_pos[dup_full_df$peptide_id == x[2]],
      strandB= dup_full_df$strand[dup_full_df$peptide_id == x[2]],
      stringsAsFactors = FALSE
    )
  }, simplify = FALSE)
})
tag_pairs_df <- do.call(rbind, unlist(tag_pairs_list, recursive = FALSE))
# -- ordering them by having startA < startB
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
         strandB = strandB_final)
# -- making sure all startA < startB
all(tag_pairs_df$startA < tag_pairs_df$startB) # -- there is only one that is =

# -- adding orientation type of each pair
tag_pairs_df <- tag_pairs_df %>%
  rowwise() %>%
  mutate(orientation = case_when(
    strandA == strandB ~ "tandem",
    strandA == 1 & strandB == -1 ~ "convergent",
    strandA == -1 & strandB == 1 ~ "divergent",
    TRUE ~ NA_character_  # in case there are unexpected strand values
  )) %>%
  ungroup()

write.table(tag_pairs_df, '../../output/duplication_classes/TAGs/TAG_gene_pairs.tsv', sep='\t', row.names=FALSE, quote=FALSE)


# -- orientation of TAG pairs with %
orientation_counts <- tag_pairs_df %>%
  group_by(orientation) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(percentage = (count / sum(count)) * 100)
print(orientation_counts)

observed <- table(tag_pairs_df$orientation)

# Suppose genome-wide proportions are known:
# tandem = 0.5, convergent = 0.25, divergent = 0.25
expected_prop <- c(tandem = 0.5, convergent = 0.25, divergent = 0.25)
expected <- sum(observed) * expected_prop[names(observed)]

chisq.test(x = observed, p = expected_prop[names(observed)])
