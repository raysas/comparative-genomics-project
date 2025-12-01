longest_family_members
# -- save it in text file


# -- distribution of family members across chromosomes
longest_family_dup_df<-dup_full_df[dup_full_df$peptide_id %in% longest_family_members, ]
table(longest_family_dup_df$chromosome)

# -- how much of them are TAGs 
longest_family_tags_df<-tags_df[tags_df$peptide_id %in% longest_family_members, ]
table(longest_family_tags_df$TAG)
length(table(longest_family_tags_df$TAG))
order(table(longest_family_tags_df$TAG), decreasing=TRUE)

# -- most of them are not (713 out of total)
# -- the rest are split to 168 different TAG blocks