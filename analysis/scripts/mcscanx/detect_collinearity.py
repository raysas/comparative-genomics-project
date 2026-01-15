#!/usr/bin/env python3
"""
Collinearity detection and anchor extraction from GFF + Ks results.
Reimplements MCScanX anchor logic without using MCScanX.
Validates against actual MCScanX output.

Inputs:
  --gff      Path to genome GFF file
  --ks       Path to Ks results (with gene1, gene2, ks columns)
  --out      Output prefix
  --gap      Max gene gap to consider collinear (default: 10)
  --mcscanx  (Optional) MCScanX collinearity file for validation
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
import re

parser = argparse.ArgumentParser()
parser.add_argument('--gff', required=True, help='Genome GFF file')
parser.add_argument('--ks', required=True, help='Ks results TSV')
parser.add_argument('--out', required=True, help='Output prefix')
parser.add_argument('--gap', type=int, default=10, help='Max gene gap for collinearity')
parser.add_argument('--mcscanx', default=None, help='(Optional) MCScanX collinearity file for validation')
args = parser.parse_args()

out_dir = Path(args.out).parent
out_dir.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("COLLINEARITY DETECTION & ANCHOR EXTRACTION (Non-MCScanX Implementation)")
print("=" * 80)

# ============================================================================
# STEP 1: Parse GFF to get gene order, chromosome, position, orientation
# ============================================================================
print(f"\n[1/5] Parsing GFF: {args.gff}")

genes = {}  # transcript_id (KRH*) -> {chrom, start, end, strand}
chrom_genes = defaultdict(list)  # chrom -> [(transcript_id, start, end, strand), ...]
transcript_bounds = defaultdict(lambda: {'min_start': float('inf'), 'max_end': 0})

gff_lines_read = 0
with open(args.gff, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        
        gff_lines_read += 1
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        
        chrom = parts[0]
        ftype = parts[2]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attr = parts[8]
        
        # Look for any KRH ID in attributes (simple regex)
        krh_match = re.search(r'(KRH\d+)', attr)
        if not krh_match:
            continue
        
        transcript_id = krh_match.group(1)
        
        # Update transcript bounds
        if transcript_id not in transcript_bounds[transcript_id]:
            transcript_bounds[transcript_id]['chrom'] = chrom
            transcript_bounds[transcript_id]['strand'] = strand
        
        transcript_bounds[transcript_id]['min_start'] = min(transcript_bounds[transcript_id]['min_start'], start)
        transcript_bounds[transcript_id]['max_end'] = max(transcript_bounds[transcript_id]['max_end'], end)

print(f"  GFF lines read: {gff_lines_read}")

# Finalize gene coordinates
for transcript_id, bounds in transcript_bounds.items():
    if 'chrom' in bounds:
        genes[transcript_id] = {
            'chrom': bounds['chrom'],
            'start': bounds['min_start'],
            'end': bounds['max_end'],
            'strand': bounds['strand']
        }
        chrom_genes[bounds['chrom']].append((transcript_id, bounds['min_start'], bounds['max_end'], bounds['strand']))

# Sort genes by position on each chromosome
for chrom in chrom_genes:
    chrom_genes[chrom].sort(key=lambda x: x[1])

print(f"  Loaded {len(genes)} transcripts from {len(chrom_genes)} chromosomes")

# Sort genes by position on each chromosome
for chrom in chrom_genes:
    chrom_genes[chrom].sort(key=lambda x: x[1])
    # Add order index
    for idx, (gene_id, start, end, strand) in enumerate(chrom_genes[chrom]):
        genes[gene_id]['order'] = idx

# ============================================================================
# STEP 2: Load Ks results and build homology map
# ============================================================================
print(f"\n[2/5] Loading Ks results: {args.ks}")

ks_df = pd.read_csv(args.ks, sep='\t', usecols=['gene1', 'gene2', 'ks'])
print(f"  Loaded {len(ks_df)} gene pairs")

# Build bidirectional homology map: gene -> set of homologs (for fast lookup)
homology = defaultdict(set)
ks_map = {}  # (g1, g2) -> ks (normalized key)

for _, row in ks_df.iterrows():
    g1 = row['gene1']
    g2 = row['gene2']
    ks = row['ks']
    
    if g1 in genes and g2 in genes:
        homology[g1].add(g2)
        homology[g2].add(g1)
        # Store Ks value with normalized key
        key = tuple(sorted([g1, g2]))
        ks_map[key] = ks

print(f"  Built homology map for {len(homology)} genes")

# ============================================================================
# STEP 3: Detect collinear blocks
# ============================================================================
print(f"\n[3/5] Detecting collinear blocks (gap threshold: {args.gap} genes)")

def find_collinear_blocks(chrom, genes_list, homology, gap_threshold):
    """Find collinear blocks on a single chromosome."""
    blocks = []
    visited = set()
    
    for i, (gene, start, end, strand) in enumerate(genes_list):
        if gene in visited or gene not in homology:
            continue
        
        # Try to extend from this gene
        block_positions = [i]
        block_genes = [gene]
        visited.add(gene)
        
        # Look for collinear chain: find homologs nearby that form a chain
        j = i + 1
        while j < len(genes_list) and j - i <= gap_threshold * 2:
            next_gene, _, _, next_strand = genes_list[j]
            
            # Check if next_gene is homologous to any recent gene in chain
            for chain_gene in block_genes[-gap_threshold:]:
                if next_gene in homology.get(chain_gene, set()):
                    # Accept if orientation is consistent or within tandem tolerance
                    block_positions.append(j)
                    block_genes.append(next_gene)
                    visited.add(next_gene)
                    break
            
            j += 1
        
        if len(block_genes) >= 2:  # Only keep blocks with 2+ pairs
            blocks.append((block_genes, block_positions))
    
    return blocks

collinear_blocks = defaultdict(list)
for chrom, genes_list in chrom_genes.items():
    blocks = find_collinear_blocks(chrom, genes_list, homology, args.gap)
    collinear_blocks[chrom] = blocks

total_blocks = sum(len(b) for b in collinear_blocks.values())
total_genes_in_blocks = sum(sum(len(genes_list) for genes_list, _ in blocks) for blocks in collinear_blocks.values())
print(f"  Found {total_blocks} collinear blocks with {total_genes_in_blocks} genes")

# ============================================================================
# STEP 4: Extract anchor pairs from collinear blocks
# ============================================================================
print(f"\n[4/5] Extracting anchor pairs from collinear blocks")

anchor_pairs = []

for chrom, blocks in collinear_blocks.items():
    for block_genes, block_positions in blocks:
        # For each gene in the block, find pairs with other block genes
        for i in range(len(block_genes)):
            for j in range(i+1, min(i + args.gap + 1, len(block_genes))):  # Only nearby genes
                g1 = block_genes[i]
                g2 = block_genes[j]
                
                # Get Ks value
                key = tuple(sorted([g1, g2]))
                ks_val = ks_map.get(key)
                
                if ks_val is not None:
                    anchor_pairs.append({
                        'gene1': g1,
                        'gene2': g2,
                        'ks': ks_val,
                        'chrom': chrom,
                        'collinear': True
                    })

anchor_df = pd.DataFrame(anchor_pairs)
print(f"  Extracted {len(anchor_df)} anchor pairs")

# Save anchor pairs
anchor_file = f"{args.out}_anchors.tsv"
anchor_df.to_csv(anchor_file, sep='\t', index=False)
print(f"  Saved to: {anchor_file}")

# ============================================================================
# STEP 5: Validation against MCScanX (if provided)
# ============================================================================
if args.mcscanx:
    print(f"\n[5/5] Validating against MCScanX: {args.mcscanx}")
    
    # Parse MCScanX collinearity file
    mcscanx_pairs = set()
    mcscanx_blocks = 0
    
    with open(args.mcscanx, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            
            # MCScanX format: "  0-  0:	KRH74855	KRH76270	  3e+01"
            # Skip lines without colons (header/stats lines)
            if ':' not in line:
                continue
            
            # Split by tabs after the colon part
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            # After the index (e.g., "0-  0:"), we have gene1 gene2 and optionally score
            # parts[0] might be "0-  0:" or similar
            # parts[1] should be gene1
            # parts[2] should be gene2
            
            try:
                g1 = parts[1].strip()
                g2 = parts[2].strip()
                
                # Validate that these are gene IDs (contain KRH)
                if 'KRH' not in g1 or 'KRH' not in g2:
                    continue
                
                # Normalize pair order
                pair = tuple(sorted([g1, g2]))
                mcscanx_pairs.add(pair)
            except (IndexError, ValueError):
                continue
    
    # Compare
    our_pairs = set(tuple(sorted([row['gene1'], row['gene2']])) for _, row in anchor_df.iterrows())
    
    overlap = our_pairs & mcscanx_pairs
    mcscanx_only = mcscanx_pairs - our_pairs
    ours_only = our_pairs - mcscanx_pairs
    
    print(f"\n  MCScanX anchors: {len(mcscanx_pairs)}")
    print(f"  Our anchors:    {len(our_pairs)}")
    print(f"  Overlap:        {len(overlap)} ({100*len(overlap)/max(len(mcscanx_pairs), 1):.1f}% of MCScanX)")
    print(f"  MCScanX-only:   {len(mcscanx_only)}")
    print(f"  Ours-only:      {len(ours_only)}")
    
    # Save comparison
    validation_file = f"{args.out}_validation.txt"
    with open(validation_file, 'w') as f:
        f.write("COLLINEARITY DETECTION VALIDATION\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"MCScanX anchors: {len(mcscanx_pairs)}\n")
        f.write(f"Our anchors:    {len(our_pairs)}\n")
        f.write(f"Overlap:        {len(overlap)} ({100*len(overlap)/max(len(mcscanx_pairs), 1):.1f}%)\n")
        f.write(f"MCScanX-only:   {len(mcscanx_only)}\n")
        f.write(f"Ours-only:      {len(ours_only)}\n")
    
    print(f"\n  Validation results saved to: {validation_file}")

print("\n" + "=" * 80)
print("COMPLETE")
print("=" * 80)
