#!/usr/bin/env python3
"""
Identify Tandemly Arrayed Genes (TAGs)
"""

import pandas as pd
import numpy as np
import argparse
from collections import defaultdict

def load_gene_positions(gff_file):
    """Load gene positions from GFF file"""
    genes = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            
            # Filter for genes only
            if feature_type != 'gene':
                continue
            
            # Parse gene ID from attributes
            attrs = {}
            for item in attributes.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    attrs[key.strip()] = value.strip()
            
            gene_id = attrs.get('ID', attrs.get('gene_id', ''))
            
            genes.append({
                'Gene_ID': gene_id,
                'Chromosome': seqid,
                'Start': int(start),
                'End': int(end),
                'Strand': strand
            })
    
    return pd.DataFrame(genes)

def identify_tags_distance_based(gene_families, gene_positions, max_distance_bp=100000):
    """
    Identify TAGs based on physical distance
    Two genes in same family are TAGs if on same chromosome and within max_distance_bp
    """
    
    print(f"Identifying TAGs with distance criterion: {max_distance_bp} bp...")
    
    # Merge gene families with positions
    df = gene_families.merge(gene_positions, on='Gene_ID', how='inner')
    
    # Sort by chromosome and position
    df = df.sort_values(['Chromosome', 'Start'])
    
    tag_pairs = []
    tag_clusters = defaultdict(set)
    
    # Group by family
    for family_id, family_df in df.groupby('Family_ID'):
        family_genes = family_df.to_dict('records')
        
        # Compare all pairs within family
        for i in range(len(family_genes)):
            for j in range(i+1, len(family_genes)):
                gene1 = family_genes[i]
                gene2 = family_genes[j]
                
                # Check if on same chromosome
                if gene1['Chromosome'] != gene2['Chromosome']:
                    continue
                
                # Calculate distance (between closest ends)
                distance = abs(gene1['Start'] - gene2['Start'])
                
                # Check if within distance threshold
                if distance <= max_distance_bp:
                    tag_pairs.append({
                        'Gene1': gene1['Gene_ID'],
                        'Gene2': gene2['Gene_ID'],
                        'Family_ID': family_id,
                        'Chromosome': gene1['Chromosome'],
                        'Distance_bp': distance,
                        'Gene1_Start': gene1['Start'],
                        'Gene2_Start': gene2['Start'],
                        'Strand1': gene1['Strand'],
                        'Strand2': gene2['Strand'],
                        'Same_Orientation': gene1['Strand'] == gene2['Strand']
                    })
                    
                    # Add to cluster
                    tag_clusters[family_id].add(gene1['Gene_ID'])
                    tag_clusters[family_id].add(gene2['Gene_ID'])
    
    return pd.DataFrame(tag_pairs), tag_clusters

def identify_tags_gene_count(gene_families, gene_positions, max_genes_apart=10):
    """
    Identify TAGs based on gene count
    Two genes in same family are TAGs if within max_genes_apart on same chromosome
    """
    
    print(f"Identifying TAGs with gene count criterion: {max_genes_apart} genes apart...")
    
    # Merge and sort
    df = gene_families.merge(gene_positions, on='Gene_ID', how='inner')
    df = df.sort_values(['Chromosome', 'Start'])
    
    # Add gene rank within chromosome
    df['Gene_Rank'] = df.groupby('Chromosome').cumcount() + 1
    
    tag_pairs = []
    tag_clusters = defaultdict(set)
    
    # Group by family
    for family_id, family_df in df.groupby('Family_ID'):
        family_genes = family_df.to_dict('records')
        
        # Compare all pairs
        for i in range(len(family_genes)):
            for j in range(i+1, len(family_genes)):
                gene1 = family_genes[i]
                gene2 = family_genes[j]
                
                # Check if on same chromosome
                if gene1['Chromosome'] != gene2['Chromosome']:
                    continue
                
                # Calculate gene distance
                gene_distance = abs(gene1['Gene_Rank'] - gene2['Gene_Rank'])
                
                # Check if within gene count threshold
                if gene_distance <= max_genes_apart:
                    tag_pairs.append({
                        'Gene1': gene1['Gene_ID'],
                        'Gene2': gene2['Gene_ID'],
                        'Family_ID': family_id,
                        'Chromosome': gene1['Chromosome'],
                        'Genes_Apart': gene_distance,
                        'Distance_bp': abs(gene1['Start'] - gene2['Start']),
                        'Strand1': gene1['Strand'],
                        'Strand2': gene2['Strand'],
                        'Same_Orientation': gene1['Strand'] == gene2['Strand']
                    })
                    
                    tag_clusters[family_id].add(gene1['Gene_ID'])
                    tag_clusters[family_id].add(gene2['Gene_ID'])
    
    return pd.DataFrame(tag_pairs), tag_clusters

def main():
    parser = argparse.ArgumentParser(description='Identify Tandemly Arrayed Genes (TAGs)')
    parser.add_argument('--families', default='../../output/clusters/protein_families_network_50families_max5.tsv',
                       help='Gene families file (TSV: Gene_ID, Family_ID)')
    parser.add_argument('--gff', default='../../data/annotations/Glycine_max.gff3',
                       help='Gene annotation GFF file')
    parser.add_argument('--gene-pos', default='../../data/reference/gene_positions.tsv',
                       help='Pre-computed gene positions file (optional, will parse GFF if not provided)')
    parser.add_argument('--method', choices=['distance', 'gene_count', 'both'], default='both',
                       help='TAG identification method')
    parser.add_argument('--max-distance', type=int, default=100000,
                       help='Maximum distance in bp for distance-based method (default: 100kb)')
    parser.add_argument('--max-genes', type=int, default=10,
                       help='Maximum genes apart for gene-count method (default: 10)')
    parser.add_argument('--output-dir', default='../../output/statistics/',
                       help='Output directory')
    
    args = parser.parse_args()
    
    # Load gene families
    print(f"Loading gene families from {args.families}...")
    gene_families = pd.read_csv(args.families, sep='\t', names=['Gene_ID', 'Family_ID'])
    print(f"  Loaded {len(gene_families)} genes in {gene_families['Family_ID'].nunique()} families")
    
    # Load gene positions
    if args.gene_pos and pd.io.common.file_exists(args.gene_pos):
        print(f"Loading gene positions from {args.gene_pos}...")
        gene_positions = pd.read_csv(args.gene_pos, sep='\t')
    else:
        print(f"Parsing gene positions from GFF: {args.gff}...")
        gene_positions = load_gene_positions(args.gff)
        # Save for future use
        gene_positions.to_csv(f"{args.output_dir}/gene_positions.tsv", sep='\t', index=False)
    
    print(f"  Loaded positions for {len(gene_positions)} genes")
    
    # Identify TAGs
    if args.method in ['distance', 'both']:
        print("\n=== Distance-based TAG identification ===")
        tag_pairs_dist, tag_clusters_dist = identify_tags_distance_based(
            gene_families, gene_positions, args.max_distance
        )
        
        # Save results
        output_file = f"{args.output_dir}/TAGs_distance_{args.max_distance}bp.tsv"
        tag_pairs_dist.to_csv(output_file, sep='\t', index=False)
        print(f"  Identified {len(tag_pairs_dist)} TAG pairs")
        print(f"  Total TAG genes: {sum(len(cluster) for cluster in tag_clusters_dist.values())}")
        print(f"  Saved to: {output_file}")
        
        # Orientation statistics
        same_orient = tag_pairs_dist['Same_Orientation'].sum()
        total = len(tag_pairs_dist)
        print(f"  Same orientation: {same_orient}/{total} ({same_orient/total*100:.1f}%)")
    
    if args.method in ['gene_count', 'both']:
        print("\n=== Gene-count-based TAG identification ===")
        tag_pairs_genes, tag_clusters_genes = identify_tags_gene_count(
            gene_families, gene_positions, args.max_genes
        )
        
        # Save results
        output_file = f"{args.output_dir}/TAGs_genecount_{args.max_genes}genes.tsv"
        tag_pairs_genes.to_csv(output_file, sep='\t', index=False)
        print(f"  Identified {len(tag_pairs_genes)} TAG pairs")
        print(f"  Total TAG genes: {sum(len(cluster) for cluster in tag_clusters_genes.values())}")
        print(f"  Saved to: {output_file}")
        
        # Orientation statistics
        same_orient = tag_pairs_genes['Same_Orientation'].sum()
        total = len(tag_pairs_genes)
        print(f"  Same orientation: {same_orient}/{total} ({same_orient/total*100:.1f}%)")
    
    # Create gene annotation file (TAG vs non-TAG)
    all_genes = gene_families['Gene_ID'].unique()
    if args.method == 'distance':
        tag_genes = set().union(*tag_clusters_dist.values())
    elif args.method == 'gene_count':
        tag_genes = set().union(*tag_clusters_genes.values())
    else:  # both - use union
        tag_genes = set().union(*tag_clusters_dist.values(), *tag_clusters_genes.values())
    
    gene_annotation = pd.DataFrame({
        'Gene_ID': all_genes,
        'Is_TAG': [gene in tag_genes for gene in all_genes]
    })
    
    output_file = f"{args.output_dir}/gene_TAG_annotation.tsv"
    gene_annotation.to_csv(output_file, sep='\t', index=False)
    print(f"\nGene TAG annotation saved to: {output_file}")
    print(f"  TAGs: {gene_annotation['Is_TAG'].sum()}")
    print(f"  Non-TAGs: {(~gene_annotation['Is_TAG']).sum()}")
    
    print("\nDone!")

if __name__ == '__main__':
    main()
