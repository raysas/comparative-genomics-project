# 1. Load file (correct your path)
# lines <- readLines("C:/Users/eyabt/Desktop/comparative-genomics-project/examples/clustering_output.txt")
lines <- readLines("../../output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.txt")

# 2. Remove empty lines (if any)
lines <- lines[nzchar(lines)]

# 3. Count number of genes per family (per line)
# Split on tabs and spaces
family_sizes <- sapply(strsplit(lines, "[\t ]+"), length)

# 4. Distribution: how many families of size 1, 2, 3...
size_distribution <- table(family_sizes)

# 5. Plot result
barplot(size_distribution,
        xlab = "Family size (number of genes)",
        ylab = "Number of families",
        main = "Gene family size distribution",
        las = 1,
        col='dodgerblue2')

