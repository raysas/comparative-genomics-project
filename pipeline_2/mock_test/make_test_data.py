#!/usr/bin/env python3
"""
Generate test dataset for Pipeline 2.
Creates balanced sample from full clustering results.
"""

import pandas as pd
import sys
import os
from pathlib import Path

def create_test_dataset(input_file, output_file, size="medium"):
    """
    Create test dataset with different size options.
    
    Args:
        input_file: Full clustering results file
        output_file: Output file for test dataset
        size: "small" (~50 genes), "medium" (~150 genes), "large" (~300 genes)
    """
    
    try:
        df = pd.read_csv(input_file, sep="\t")
    except FileNotFoundError:
        print(f"❌ Input file not found: {input_file}")
        return False
    
    # Get family size distribution
    family_sizes = df.groupby('family').size()
    
    # Define selection criteria based on size
    if size == "small":
        # Small test: 3-4 families, ~50 genes total
        small_fams = family_sizes[family_sizes.between(2, 4)].index[:2]
        medium_fams = family_sizes[family_sizes.between(5, 8)].index[:2]
        large_fams = []
        print("Creating SMALL test dataset (~50 genes)")
        
    elif size == "large":
        # Large test: 10+ families, ~300 genes total
        small_fams = family_sizes[family_sizes.between(2, 4)].index[:8]
        medium_fams = family_sizes[family_sizes.between(5, 10)].index[:6]
        large_fams = family_sizes[family_sizes.between(10, 20)].index[:3]
        print("Creating LARGE test dataset (~300 genes)")
        
    else:  # medium (default)
        # Medium test: 6-8 families, ~150 genes total
        small_fams = family_sizes[family_sizes.between(2, 4)].index[:3]
        medium_fams = family_sizes[family_sizes.between(5, 10)].index[:3]
        large_fams = family_sizes[family_sizes.between(10, 20)].index[:2]
        print("Creating MEDIUM test dataset (~150 genes)")
    
    # Combine selected families
    selected_families = list(small_fams) + list(medium_fams) + list(large_fams)
    
    if not selected_families:
        print("❌ No suitable families found")
        return False
    
    # Extract and renumber families
    test_data = df[df['family'].isin(selected_families)].copy()
    family_mapping = {old: new for new, old in enumerate(selected_families, 1)}
    test_data['family'] = test_data['family'].map(family_mapping)
    
    # Save dataset
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    test_data.to_csv(output_file, sep="\t", index=False)
    
    # Print summary
    total_pairs = sum(len(group) * (len(group) - 1) // 2 
                     for _, group in test_data.groupby('family'))
    
    print(f"✅ Created: {len(test_data)} genes, {len(selected_families)} families, {total_pairs} pairs")
    return True

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: make_test_data.py <input_file> <output_file> [small|medium|large]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] 
    size = sys.argv[3] if len(sys.argv) > 3 else "medium"
    
    success = create_test_dataset(input_file, output_file, size)
    sys.exit(0 if success else 1)