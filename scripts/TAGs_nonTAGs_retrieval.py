#!/bin/env python3

stringency='high'

total_genes_file=f'output/gene_lists/duplicated_genes/duplicated_genes_{stringency}.txt'
tag_genes_file=f'output/gene_lists/TAGs/spacer_based/TAGs_{stringency}.txt'

with open(total_genes_file, 'r') as f:
    total_genes = set(line.strip() for line in f)

with open(tag_genes_file, 'r') as f:
    tag_genes = set(line.strip() for line in f)

non_tag_genes = total_genes - tag_genes
output_file = f'output/gene_lists/nonTAGs/non_TAGs_{stringency}.txt'

with open(output_file, 'w') as f:
    for gene in sorted(non_tag_genes):
        f.write(gene + '\n')