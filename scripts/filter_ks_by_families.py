#!/usr/bin/env python3
"""
Filter Ks results based on gene families from filtered BLAST results.
Only keeps gene pairs where both genes are present in the filtered family list
and belong to the same family.

Usage:
    python3 filter_ks_by_families.py ks_results.tsv family_filter.tsv output.tsv
    
Example:
    python3 scripts/filter_ks_by_families.py \
        output/ks_results/ks_results_filtered.tsv \
        output/filtered/protein_families_filtered_blast_results_id50_qcov70_scov70_wcol12_network.tsv \
        output/ks_results/ks_results_id50_qcov70_scov70.tsv
"""

import pandas as pd
import sys
import os

def filter_ks_by_families(ks_file, family_file, output_file):
    """
    Filter Ks results to include only gene pairs where both genes
    are in the filtered family list and belong to the same family.
    
    Args:
        ks_file: Path to Ks results file (gene1, gene2, family, ks, ka, etc.)
        family_file: Path to filtered family assignments (geneName, family)
        output_file: Path to output filtered Ks results
    """
    
    # Load data
    print(f"Loading Ks results: {ks_file}")
    ks_df = pd.read_csv(ks_file, sep='\t')
    print(f"  Total Ks pairs: {len(ks_df)}")
    
    print(f"\nLoading family filter: {family_file}")
    family_df = pd.read_csv(family_file, sep='\t')
    print(f"  Total genes in filter: {len(family_df)}")
    print(f"  Total families: {family_df['family'].nunique()}")
    
    # Create gene to family mapping
    gene_to_family = dict(zip(family_df['geneName'], family_df['family']))
    
    # Filter Ks results
    print("\nFiltering Ks results...")
    
    def is_valid_pair(row):
        gene1 = row['gene1']
        gene2 = row['gene2']
        
        # Check if both genes are in the filtered list
        if gene1 not in gene_to_family or gene2 not in gene_to_family:
            return False
        
        # Check if they belong to the same family
        if gene_to_family[gene1] != gene_to_family[gene2]:
            return False
        
        return True
    
    ks_filtered = ks_df[ks_df.apply(is_valid_pair, axis=1)].copy()
    
    print(f"  Filtered Ks pairs: {len(ks_filtered)}")
    print(f"  Retained: {len(ks_filtered) / len(ks_df) * 100:.2f}%")
    
    # Save results
    print(f"\nSaving filtered results to: {output_file}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    ks_filtered.to_csv(output_file, sep='\t', index=False)
    
    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Original pairs: {len(ks_df)}")
    print(f"Filtered pairs: {len(ks_filtered)}")
    print(f"Removed pairs: {len(ks_df) - len(ks_filtered)}")
    print(f"Families represented: {ks_filtered['family'].nunique()}")
    
    if 'ks' in ks_filtered.columns:
        print(f"\nKs statistics (filtered):")
        print(f"  Mean: {ks_filtered['ks'].mean():.4f}")
        print(f"  Median: {ks_filtered['ks'].median():.4f}")
        print(f"  Min: {ks_filtered['ks'].min():.4f}")
        print(f"  Max: {ks_filtered['ks'].max():.4f}")
    
    return ks_filtered

def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    
    ks_file = sys.argv[1]
    family_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Check input files exist
    if not os.path.exists(ks_file):
        print(f"ERROR: Ks results file not found: {ks_file}")
        sys.exit(1)
    
    if not os.path.exists(family_file):
        print(f"ERROR: Family filter file not found: {family_file}")
        sys.exit(1)
    
    # Run filtering
    filter_ks_by_families(ks_file, family_file, output_file)
    
    print("\n✓ Filtering complete!")

if __name__ == "__main__":
    main()
