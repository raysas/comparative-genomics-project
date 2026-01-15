#!/usr/bin/env python3
"""
Infer centromere locations from gene distribution.

Centromeres appear as gene deserts (regions with few/no genes).
This script identifies the largest gap in gene distribution on each chromosome
and marks it as the likely centromere region.

Input: GFF file or protein info CSV
Output: TSV file with inferred centromere locations

Method:
  1. For each chromosome, get all gene positions
  2. Sort genes by position
  3. Find the largest gap between consecutive genes
  4. Mark that gap as the centromere (with 4Mb buffer on each side for safety)
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--gff', default=None, help='GFF file')
parser.add_argument('--protein-info', default=None, help='Protein info CSV')
parser.add_argument('--output', required=True, help='Output centromere TSV')
parser.add_argument('--buffer', type=int, default=2_000_000, 
                    help='Buffer around detected gap (default 2Mb)')
args = parser.parse_args()

print("=" * 80)
print("CENTROMERE INFERENCE FROM GENE DISTRIBUTION")
print("=" * 80)

# Load gene positions
genes_by_chr = defaultdict(list)

if args.gff:
    print(f"\nLoading genes from GFF: {args.gff}")
    with open(args.gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            
            genes_by_chr[chrom].append((start, end))

elif args.protein_info:
    print(f"\nLoading genes from protein info: {args.protein_info}")
    df = pd.read_csv(args.protein_info)
    
    for _, row in df.iterrows():
        chrom = row['chromosome']
        start = row['start_pos']
        end = row['end_pos']
        genes_by_chr[chrom].append((start, end))

else:
    print("ERROR: Must provide --gff or --protein-info")
    exit(1)

print(f"  ✓ Loaded {sum(len(genes) for genes in genes_by_chr.values())} genes")
print(f"  ✓ Across {len(genes_by_chr)} chromosomes")

# Infer centromeres
print("\nInferring centromere locations...\n")

centromeres = []
chr_summaries = []

for chrom in sorted(genes_by_chr.keys()):
    genes = genes_by_chr[chrom]
    
    if len(genes) < 3:
        print(f"  ⚠ {chrom}: Too few genes ({len(genes)}) to infer centromere")
        continue
    
    # Get midpoint of each gene
    gene_positions = [(start + end) / 2 for start, end in genes]
    gene_positions.sort()
    
    # Find largest gap
    gaps = []
    for i in range(len(gene_positions) - 1):
        gap_size = gene_positions[i + 1] - gene_positions[i]
        gap_start = gene_positions[i]
        gap_end = gene_positions[i + 1]
        gaps.append({
            'size': gap_size,
            'start': gap_start,
            'end': gap_end,
            'midpoint': (gap_start + gap_end) / 2
        })
    
    largest_gap = max(gaps, key=lambda x: x['size'])
    
    # Mark centromere with buffer
    buffer = args.buffer
    centro_start = int(max(0, largest_gap['start'] - buffer))
    centro_end = int(largest_gap['end'] + buffer)
    gap_size_mb = largest_gap['size'] / 1_000_000
    
    # Get chromosome length
    chr_length = max(pos for pos in gene_positions)
    
    centromeres.append({
        'chromosome': chrom,
        'start': centro_start,
        'end': centro_end,
        'gap_size_bp': int(largest_gap['size']),
        'gap_center_bp': int(largest_gap['midpoint']),
        'buffer_bp': buffer
    })
    
    chr_summaries.append({
        'Chromosome': chrom,
        'Genes': len(genes),
        'Largest_Gap_Mb': f"{gap_size_mb:.2f}",
        'Centromere_Start_Mb': f"{centro_start/1_000_000:.2f}",
        'Centromere_End_Mb': f"{centro_end/1_000_000:.2f}",
        'Chr_Length_Mb': f"{chr_length/1_000_000:.2f}",
        'Centromere_Position_%': f"{(largest_gap['midpoint']/chr_length)*100:.1f}"
    })
    
    print(f"  {chrom}:")
    print(f"    Genes: {len(genes)}")
    print(f"    Largest gap: {gap_size_mb:.2f} Mb")
    print(f"    Inferred centromere: {centro_start:,} - {centro_end:,} bp")
    print(f"    Position on chr: {(largest_gap['midpoint']/chr_length)*100:.1f}%")

# Save results
print(f"\nSaving centromere predictions...")

# Save TSV format
df_centromeres = pd.DataFrame(centromeres)
df_centromeres.to_csv(args.output, sep='\t', index=False)
print(f"  ✓ Saved to: {args.output}")

# Save summary
summary_path = Path(args.output).parent / f"{Path(args.output).stem}_summary.txt"
with open(summary_path, 'w') as f:
    f.write("INFERRED CENTROMERE LOCATIONS\n")
    f.write("=" * 100 + "\n")
    f.write(f"Method: Largest gap in gene distribution (with {args.buffer/1_000_000:.1f}Mb buffer)\n\n")
    
    summary_df = pd.DataFrame(chr_summaries)
    f.write(summary_df.to_string(index=False))
    f.write("\n\n" + "=" * 100 + "\n")
    f.write("\nNotes:\n")
    f.write("  - Centromere position should coincide with the largest gap in gene distribution\n")
    f.write("  - Buffer extends predictions by specified amount on both sides\n")
    f.write("  - Verify with known centromere positions for your organism if available\n")
    f.write("  - These are rough estimates; actual centromeres may vary\n")

print(f"  ✓ Saved summary to: {summary_path}")

# Print summary table
print("\nSummary:")
print(pd.DataFrame(chr_summaries).to_string(index=False))

print("\n" + "=" * 80)
print("✓ COMPLETE")
print("=" * 80)
