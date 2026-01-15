#!/usr/bin/env python3
"""
Analyze correlation between WGD/TAG duplicates and chromosomal structural patterns.

Investigates:
  1. Proximity to chromosome ends (telomeres)
  2. Clustering patterns (hotspots of duplication)
  3. Correlation with gene density
  4. Pericentromeric regions (if data available)
  5. Distribution relative to chromosome length
  6. Co-localization of different duplication types

Usage:
    python3 analyze_structural_patterns.py \
        --classified gene_pairs_classified.tsv \
        --protein-info protein_info_longest.csv \
        --outdir analysis/structural_analysis
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import chi2_contingency, kstest, spearmanr
from collections import defaultdict

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150

# Define structural regions
TELOMERE_DISTANCE = 4_000_000  # 4Mb from chromosome ends
PERICENTROMERIC_DISTANCE = 2_000_000  # 2Mb from centromere (conservative)
CLUSTERING_WINDOW = 500_000  # 500kb window for clustering analysis


def load_data(classified_path: Path, protein_info_path: Path, centromere_path: Path = None):
    """Load classified gene pairs, protein metadata, and optionally centromere locations."""
    print("Loading data...")
    
    classified_df = pd.read_csv(classified_path, sep='\t')
    protein_df = pd.read_csv(protein_info_path)
    
    print(f"  ✓ Loaded {len(classified_df)} classified gene pairs")
    print(f"  ✓ Loaded {len(protein_df)} genes")
    
    centromere_df = None
    if centromere_path and Path(centromere_path).exists():
        centromere_df = pd.read_csv(centromere_path, sep='\t')
        print(f"  ✓ Loaded {len(centromere_df)} centromere regions")
    
    return classified_df, protein_df, centromere_df


def get_chromosome_lengths(protein_df):
    """Calculate chromosome lengths from gene positions."""
    print("Calculating chromosome lengths...")
    
    chr_lengths = {}
    for chrom in protein_df['chromosome'].unique():
        if pd.isna(chrom):
            continue
        # Convert to string for consistent key format
        chrom_str = str(int(chrom))
        chr_genes = protein_df[protein_df['chromosome'] == chrom]
        max_pos = chr_genes['end_pos'].max()
        chr_lengths[chrom_str] = max_pos
    
    print(f"  ✓ Calculated lengths for {len(chr_lengths)} chromosomes")
    return chr_lengths


def classify_telomeric_region(pos, chr_length, distance=TELOMERE_DISTANCE):
    """Check if position is in telomeric region (4Mb from ends)."""
    if pd.isna(pos) or pd.isna(chr_length):
        return 'Unknown'
    
    if pos <= distance:
        return 'Telomere-proximal'
    elif pos >= chr_length - distance:
        return 'Telomere-proximal'
    else:
        return 'Internal'


def classify_pericentromeric_region(pos, centromere_start, centromere_end, distance=PERICENTROMERIC_DISTANCE):
    """Check if position is in pericentromeric region (4Mb from centromere)."""
    if pd.isna(pos) or pd.isna(centromere_start) or pd.isna(centromere_end):
        return 'Unknown'
    
    # Check distance from centromere (use closest edge)
    dist_to_center = min(abs(pos - centromere_start), abs(pos - centromere_end))
    
    if dist_to_center <= distance:
        return 'Pericentromeric'
    else:
        return 'Internal'


def classify_position_zone(pos, chr_length):
    """Classify position into quintiles of chromosome."""
    if pd.isna(pos) or pd.isna(chr_length):
        return 'Unknown'
    
    relative_pos = pos / chr_length
    
    if relative_pos < 0.2:
        return 'Q1 (0-20%)'
    elif relative_pos < 0.4:
        return 'Q2 (20-40%)'
    elif relative_pos < 0.6:
        return 'Q3 (40-60%)'
    elif relative_pos < 0.8:
        return 'Q4 (60-80%)'
    else:
        return 'Q5 (80-100%)'


def add_structural_features(classified_df, protein_df, chr_lengths, centromere_df=None):
    """Add structural feature annotations to gene pairs."""
    print("Adding structural feature annotations...")
    
    # Create centromere lookup
    centromere_map = {}
    if centromere_df is not None:
        for _, row in centromere_df.iterrows():
            chrom = row['chromosome'] if 'chromosome' in row else row.get('chrom')
            # Convert chromosome to string for consistent matching
            chrom = str(int(chrom)) if not pd.isna(chrom) else None
            if chrom:
                centromere_map[chrom] = {
                    'start': row['start'] if 'start' in row else row.get('centromere_start'),
                    'end': row['end'] if 'end' in row else row.get('centromere_end')
                }
    
    # Create gene info lookup
    gene_info = protein_df.set_index('peptide_id')[
        ['chromosome', 'start_pos', 'end_pos']
    ].to_dict('index')
    
    # Add features for each gene
    structural_features = []
    
    for _, row in classified_df.iterrows():
        gene1 = row['gene1']
        gene2 = row['gene2']
        
        features = {
            'pair_id': f"{gene1}_{gene2}",
            'duplication_type': row['duplication_type'],
            'ks': row['ks']
        }
        
        # Get gene positions
        if gene1 in gene_info and gene2 in gene_info:
            info1 = gene_info[gene1]
            info2 = gene_info[gene2]
            
            chr1 = str(int(info1['chromosome'])) if not pd.isna(info1['chromosome']) else None
            chr2 = str(int(info2['chromosome'])) if not pd.isna(info2['chromosome']) else None
            pos1 = info1['start_pos']
            pos2 = info2['start_pos']
            
            # Same chromosome info
            features['same_chromosome'] = (chr1 == chr2)
            
            # Telomere proximity for each gene
            if chr1 in chr_lengths:
                features['gene1_telomere_region'] = classify_telomeric_region(
                    pos1, chr_lengths[chr1]
                )
                features['gene1_zone'] = classify_position_zone(
                    pos1, chr_lengths[chr1]
                )
            
            if chr2 in chr_lengths:
                features['gene2_telomere_region'] = classify_telomeric_region(
                    pos2, chr_lengths[chr2]
                )
                features['gene2_zone'] = classify_position_zone(
                    pos2, chr_lengths[chr2]
                )
            
            # Pericentromeric regions (if centromere data available)
            if centromere_df is not None:
                if chr1 in centromere_map:
                    features['gene1_pericentromeric_region'] = classify_pericentromeric_region(
                        pos1, centromere_map[chr1]['start'], centromere_map[chr1]['end']
                    )
                else:
                    features['gene1_pericentromeric_region'] = 'Unknown'
                
                if chr2 in centromere_map:
                    features['gene2_pericentromeric_region'] = classify_pericentromeric_region(
                        pos2, centromere_map[chr2]['start'], centromere_map[chr2]['end']
                    )
                else:
                    features['gene2_pericentromeric_region'] = 'Unknown'
            
            # Both in telomeric regions?
            if features.get('gene1_telomere_region') == 'Telomere-proximal' and \
               features.get('gene2_telomere_region') == 'Telomere-proximal':
                features['both_telomeric'] = True
            else:
                features['both_telomeric'] = False
            
            # Both in pericentromeric regions?
            if centromere_df is not None:
                if features.get('gene1_pericentromeric_region') == 'Pericentromeric' and \
                   features.get('gene2_pericentromeric_region') == 'Pericentromeric':
                    features['both_pericentromeric'] = True
                else:
                    features['both_pericentromeric'] = False
            else:
                features['both_pericentromeric'] = False
            
            # Distance to chromosome end (for same chr pairs)
            if chr1 == chr2 and chr1 in chr_lengths:
                dist_to_start = min(pos1, pos2)
                dist_to_end = min(chr_lengths[chr1] - pos1, chr_lengths[chr1] - pos2)
                features['min_dist_to_end'] = min(dist_to_start, dist_to_end)
                features['close_to_end'] = features['min_dist_to_end'] <= TELOMERE_DISTANCE
            else:
                features['min_dist_to_end'] = np.nan
                features['close_to_end'] = False
        
        structural_features.append(features)
    
    features_df = pd.DataFrame(structural_features)
    
    print(f"  ✓ Added structural features to {len(features_df)} pairs")
    return features_df


def analyze_telomere_association(features_df, outdir: Path):
    """Analyze association between duplication types and telomeric regions."""
    print("\nAnalyzing telomere association...")
    
    # Filter to TAG and WGD only
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])].copy()
    
    # Count telomeric vs internal
    telomere_counts = pd.crosstab(
        subset['duplication_type'],
        subset['gene1_telomere_region'],
        margins=True
    )
    
    print("  Gene 1 - Telomere Association:")
    print(telomere_counts)
    
    # Chi-square test
    contingency = pd.crosstab(subset['duplication_type'], subset['gene1_telomere_region'])
    chi2, p_val, dof, expected = chi2_contingency(contingency)
    print(f"  Chi-square test: χ²={chi2:.4f}, p={p_val:.4e}")
    
    # Save results
    results = {
        'Analysis': 'Telomere Association (Gene1)',
        'Chi2': chi2,
        'P-value': p_val,
        'Significant': p_val < 0.05
    }
    
    return results, contingency


def analyze_clustering_patterns(features_df, protein_df, chr_lengths, outdir: Path):
    """Identify clustering hotspots of duplications."""
    print("\nAnalyzing clustering patterns...")
    
    # Get all duplicated genes and their positions
    gene_info = protein_df.set_index('peptide_id')[
        ['chromosome', 'start_pos']
    ].to_dict('index')
    
    all_dup_genes = set()
    for col in ['gene1', 'gene2']:
        # This would be done by reading from original classified_df
        pass
    
    clustering_analysis = defaultdict(list)
    
    # For each chromosome, identify duplication hotspots
    for chrom in protein_df['chromosome'].unique():
        # Get all duplicated genes on this chromosome
        dup_genes_chr = []
        
        for idx, row in features_df[features_df['same_chromosome']].iterrows():
            # Need access to original gene info - this is simplified
            pass
        
        if len(dup_genes_chr) > 0:
            dup_genes_chr.sort(key=lambda x: x[1])  # Sort by position
            
            # Find clusters using sliding window
            clustering_analysis[chrom] = {
                'n_dup_genes': len(dup_genes_chr),
                'span': dup_genes_chr[-1][1] - dup_genes_chr[0][1] if dup_genes_chr else 0
            }
    
    print(f"  ✓ Analyzed clustering on {len(clustering_analysis)} chromosomes")
    return clustering_analysis


def analyze_zone_distribution(features_df, outdir: Path):
    """Analyze distribution of duplications across chromosome zones."""
    print("\nAnalyzing chromosome zone distribution...")
    
    # Filter to TAG and WGD
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])].copy()
    
    # Create zone distribution table
    zone_dist = pd.crosstab(
        subset['duplication_type'],
        subset['gene1_zone'],
        margins=False
    )
    
    # Reorder columns
    zone_order = ['Q1 (0-20%)', 'Q2 (20-40%)', 'Q3 (40-60%)', 'Q4 (60-80%)', 'Q5 (80-100%)']
    zone_dist = zone_dist[[col for col in zone_order if col in zone_dist.columns]]
    
    # Normalize to percentages
    zone_dist_pct = zone_dist.div(zone_dist.sum(axis=1), axis=0) * 100
    
    print("  Zone Distribution (%):")
    print(zone_dist_pct.round(2))
    
    # Chi-square test
    chi2, p_val, dof, expected = chi2_contingency(zone_dist)
    print(f"  Chi-square test: χ²={chi2:.4f}, p={p_val:.4e}")
    
    return zone_dist_pct, chi2, p_val


def analyze_ks_by_location(features_df, outdir: Path):
    """Analyze Ks values by chromosomal location."""
    print("\nAnalyzing Ks by chromosomal location...")
    
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])].copy()
    
    # Ks by telomeric region
    telomere_ks = subset.groupby(['duplication_type', 'gene1_telomere_region'])['ks'].describe()
    print("  Ks by Telomeric Region:")
    print(telomere_ks)
    
    # Ks by chromosome zone
    zone_ks = subset.groupby(['duplication_type', 'gene1_zone'])['ks'].describe()
    print("\n  Ks by Chromosome Zone:")
    print(zone_ks)
    
    return telomere_ks, zone_ks


def plot_telomere_patterns(features_df, outdir: Path):
    """Plot association with telomeric regions - separate plots."""
    print("Generating telomere pattern plots...")
    
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])].copy()
    colors = {'TAG': '#e74c3c', 'WGD': '#2ecc71'}
    
    # 1. Stacked bar: Gene 1 telomeric status
    fig, ax = plt.subplots(figsize=(10, 6))
    telomere_data = pd.crosstab(
        subset['duplication_type'],
        subset['gene1_telomere_region'],
        normalize='index'
    ) * 100
    telomere_data.plot(kind='bar', stacked=True, ax=ax, 
                       color=['#95a5a6', '#e67e22'], alpha=0.8)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title('Gene 1: Telomeric Region Distribution', fontsize=14, fontweight='bold')
    ax.legend(title='Region', bbox_to_anchor=(1.05, 1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outdir / 'telomere_gene1_distribution.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'telomere_gene1_distribution.pdf', bbox_inches='tight')
    plt.close()
    
    # 2. Both genes in telomeric regions
    fig, ax = plt.subplots(figsize=(10, 6))
    both_telomeric = subset.groupby('duplication_type')['both_telomeric'].sum()
    both_total = subset.groupby('duplication_type').size()
    both_pct = (both_telomeric / both_total * 100).fillna(0)
    
    bars = ax.bar(both_pct.index, both_pct.values, 
                  color=[colors.get(t, 'gray') for t in both_pct.index],
                  alpha=0.8, edgecolor='black', width=0.6)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title('Both Genes in Telomeric Regions', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    for bar, count in zip(bars, both_telomeric.values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height,
               f'n={int(count)}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(outdir / 'telomere_both_genes.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'telomere_both_genes.pdf', bbox_inches='tight')
    plt.close()
    
    # 3. Box plot: Ks by telomeric region
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_data = []
    for dtype in ['TAG', 'WGD']:
        for region in ['Telomere-proximal', 'Internal']:
            subset_region = subset[(subset['duplication_type'] == dtype) & 
                                   (subset['gene1_telomere_region'] == region)]['ks'].dropna()
            if len(subset_region) > 0:
                plot_data.append({
                    'Type': dtype,
                    'Region': region,
                    'Ks': subset_region.values
                })
    
    if plot_data:
        bp_data = [d['Ks'] for d in plot_data]
        bp_labels = [f"{d['Type']}\n{d['Region']}" for d in plot_data]
        bp = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True)
        
        for i, patch in enumerate(bp['boxes']):
            dtype = plot_data[i]['Type']
            patch.set_facecolor(colors.get(dtype, 'gray'))
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Ks', fontsize=12)
        ax.set_title('Ks Distribution by Telomeric Region', fontsize=14, fontweight='bold')
        ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(outdir / 'telomere_ks_boxplot.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'telomere_ks_boxplot.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved telomere pattern plots")


def plot_zone_distribution(features_df, outdir: Path):
    """Plot distribution across chromosome zones - separate plots."""
    print("Generating chromosome zone distribution plots...")
    
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])].copy()
    colors = {'TAG': '#e74c3c', 'WGD': '#2ecc71'}
    zone_order = ['Q1 (0-20%)', 'Q2 (20-40%)', 'Q3 (40-60%)', 'Q4 (60-80%)', 'Q5 (80-100%)']
    
    # 1. Stacked bar: Absolute counts
    fig, ax = plt.subplots(figsize=(10, 6))
    zone_counts = pd.crosstab(subset['duplication_type'], subset['gene1_zone'])
    zone_counts = zone_counts[[col for col in zone_order if col in zone_counts.columns]]
    zone_counts.plot(kind='bar', stacked=False, ax=ax, width=0.8)
    ax.set_ylabel('Number of gene pairs', fontsize=12)
    ax.set_title('Distribution: Gene Pair Counts by Zone', fontsize=14, fontweight='bold')
    ax.legend(title='Zone', bbox_to_anchor=(1.05, 1), fontsize=9)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outdir / 'zone_counts.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'zone_counts.pdf', bbox_inches='tight')
    plt.close()
    
    # 2. Stacked bar: Percentages
    fig, ax = plt.subplots(figsize=(10, 6))
    zone_pct = pd.crosstab(subset['duplication_type'], subset['gene1_zone'], normalize='index') * 100
    zone_pct = zone_pct[[col for col in zone_order if col in zone_pct.columns]]
    zone_pct.plot(kind='bar', stacked=True, ax=ax, width=0.8)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title('Distribution: Gene Pair Percentages by Zone', fontsize=14, fontweight='bold')
    ax.legend(title='Zone', bbox_to_anchor=(1.05, 1), fontsize=9)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outdir / 'zone_percentages.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'zone_percentages.pdf', bbox_inches='tight')
    plt.close()
    
    # 3. Mean Ks by zone
    fig, ax = plt.subplots(figsize=(10, 6))
    zone_ks_mean = subset.groupby(['duplication_type', 'gene1_zone'])['ks'].mean().unstack()
    zone_ks_mean = zone_ks_mean[[col for col in zone_order if col in zone_ks_mean.columns]]
    zone_ks_mean.plot(ax=ax, marker='o', linewidth=2)
    ax.set_ylabel('Mean Ks', fontsize=12)
    ax.set_xlabel('Duplication Type', fontsize=12)
    ax.set_title('Mean Ks by Chromosome Zone', fontsize=14, fontweight='bold')
    ax.legend(title='Zone', bbox_to_anchor=(1.05, 1), fontsize=9)
    ax.grid(alpha=0.3)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    plt.tight_layout()
    plt.savefig(outdir / 'zone_ks_mean.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'zone_ks_mean.pdf', bbox_inches='tight')
    plt.close()
    
    # 4. Violin plot: Ks distribution across zones
    fig, ax = plt.subplots(figsize=(14, 6))
    plot_df = []
    for dtype in ['TAG', 'WGD']:
        for zone in zone_order:
            subset_zone = subset[(subset['duplication_type'] == dtype) & 
                                (subset['gene1_zone'] == zone)]['ks'].dropna()
            if len(subset_zone) > 0:
                for ks_val in subset_zone:
                    plot_df.append({'Type': dtype, 'Zone': zone, 'Ks': ks_val})
    
    if plot_df:
        plot_df = pd.DataFrame(plot_df)
        sns.violinplot(data=plot_df, x='Zone', y='Ks', hue='Type', ax=ax, 
                      palette=colors, inner='quartile')
        ax.set_xlabel('Chromosome Zone', fontsize=12)
        ax.set_ylabel('Ks', fontsize=12)
        ax.set_title('Ks Distribution across Chromosome Zones', fontsize=14, fontweight='bold')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(outdir / 'zone_ks_violin.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'zone_ks_violin.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved zone distribution plots")


def plot_pericentromeric_patterns(features_df, outdir: Path):
    """Plot association with pericentromeric regions - separate plots."""
    print("Generating pericentromeric pattern plots...")
    
    # Check if pericentromeric data exists
    if 'gene1_pericentromeric_region' not in features_df.columns:
        print("  ⚠ Centromere data not available, skipping pericentromeric analysis")
        return
    
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])].copy()
    subset = subset[subset['gene1_pericentromeric_region'] != 'Unknown']
    
    if len(subset) == 0:
        print("  ⚠ No pericentromeric data available, skipping")
        return
    
    colors = {'TAG': '#e74c3c', 'WGD': '#2ecc71'}
    
    # 1. Stacked bar: Gene 1 pericentromeric status
    fig, ax = plt.subplots(figsize=(10, 6))
    peri_data = pd.crosstab(
        subset['duplication_type'],
        subset['gene1_pericentromeric_region'],
        normalize='index'
    ) * 100
    peri_data.plot(kind='bar', stacked=True, ax=ax, 
                   color=['#95a5a6', '#e67e22'], alpha=0.8)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title('Gene 1: Pericentromeric Region Distribution (4Mb)', fontsize=14, fontweight='bold')
    ax.legend(title='Region', bbox_to_anchor=(1.05, 1))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outdir / 'pericentromeric_gene1_distribution.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'pericentromeric_gene1_distribution.pdf', bbox_inches='tight')
    plt.close()
    
    # 2. Both genes in pericentromeric regions
    fig, ax = plt.subplots(figsize=(10, 6))
    both_peri = subset.groupby('duplication_type')['both_pericentromeric'].sum()
    both_total = subset.groupby('duplication_type').size()
    both_pct = (both_peri / both_total * 100).fillna(0)
    
    bars = ax.bar(both_pct.index, both_pct.values, 
                  color=[colors.get(t, 'gray') for t in both_pct.index],
                  alpha=0.8, edgecolor='black', width=0.6)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.set_title('Both Genes in Pericentromeric Regions', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    for bar, count in zip(bars, both_peri.values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height,
               f'n={int(count)}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(outdir / 'pericentromeric_both_genes.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'pericentromeric_both_genes.pdf', bbox_inches='tight')
    plt.close()
    
    # 3. Comparison: Telomeric vs Pericentromeric
    fig, ax = plt.subplots(figsize=(10, 6))
    comparison_data = []
    for dtype in ['TAG', 'WGD']:
        subset_dtype = subset[subset['duplication_type'] == dtype]
        telomeric = (subset_dtype['gene1_telomere_region'] == 'Telomere-proximal').sum()
        pericentromeric = (subset_dtype['gene1_pericentromeric_region'] == 'Pericentromeric').sum()
        comparison_data.append({'Type': dtype, 'Telomeric': telomeric, 'Pericentromeric': pericentromeric})
    
    comp_df = pd.DataFrame(comparison_data).set_index('Type')
    comp_df.plot(kind='bar', ax=ax, alpha=0.8, 
                color=['#e74c3c', '#3498db'], edgecolor='black', width=0.6)
    ax.set_ylabel('Number of pairs', fontsize=12)
    ax.set_title('Distribution: Telomeric vs Pericentromeric Regions (4Mb)', fontsize=14, fontweight='bold')
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    ax.legend(title='Region')
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outdir / 'pericentromeric_telomeric_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'pericentromeric_telomeric_comparison.pdf', bbox_inches='tight')
    plt.close()
    
    # 4. Box plot: Ks by pericentromeric region
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_data = []
    for dtype in ['TAG', 'WGD']:
        for region in ['Pericentromeric', 'Internal']:
            subset_region = subset[(subset['duplication_type'] == dtype) & 
                                   (subset['gene1_pericentromeric_region'] == region)]['ks'].dropna()
            if len(subset_region) > 0:
                plot_data.append({
                    'Type': dtype,
                    'Region': region,
                    'Ks': subset_region.values
                })
    
    if plot_data:
        bp_data = [d['Ks'] for d in plot_data]
        bp_labels = [f"{d['Type']}\n{d['Region']}" for d in plot_data]
        bp = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True)
        
        for i, patch in enumerate(bp['boxes']):
            dtype = plot_data[i]['Type']
            patch.set_facecolor(colors.get(dtype, 'gray'))
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Ks', fontsize=12)
        ax.set_title('Ks Distribution by Pericentromeric Region', fontsize=14, fontweight='bold')
        ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(outdir / 'pericentromeric_ks_boxplot.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'pericentromeric_ks_boxplot.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved pericentromeric pattern plots")


def plot_distance_to_end_per_chromosome(features_df, classified_df, protein_df, centromere_df, outdir: Path):
    """Plot distance-to-end histograms per chromosome with centromere locations."""
    print("Generating per-chromosome distance-to-end plots...")
    
    if centromere_df is None or len(centromere_df) == 0:
        print("  ⚠ Centromere data not available, skipping per-chromosome analysis")
        return
    
    # Filter to same chromosome TAG/WGD pairs
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])].copy()
    subset = subset[subset['same_chromosome']].copy()
    
    # Map pair_id back to gene names from classified_df
    classified_df['pair_id_temp'] = classified_df['gene1'] + '_' + classified_df['gene2']
    gene_lookup = classified_df.set_index('pair_id_temp')[['gene1', 'gene2']].to_dict('index')
    
    subset['gene1'] = subset['pair_id'].map(lambda x: gene_lookup.get(x, {}).get('gene1'))
    subset = subset.dropna(subset=['gene1'])
    
    # Create gene to chromosome lookup
    gene_info = protein_df.set_index('peptide_id')['chromosome'].to_dict()
    
    # Add chromosome info
    subset['chromosome'] = subset['gene1'].map(gene_info)
    subset = subset.dropna(subset=['chromosome'])
    subset['chromosome'] = subset['chromosome'].apply(lambda x: str(int(x)))
    
    # Get chromosome lengths
    chr_lengths = {}
    for chrom in protein_df['chromosome'].unique():
        if pd.isna(chrom):
            continue
        chrom_str = str(int(chrom))
        chr_genes = protein_df[protein_df['chromosome'] == chrom]
        chr_lengths[chrom_str] = chr_genes['end_pos'].max()
    
    # Create centromere map
    centromere_map = {}
    for _, row in centromere_df.iterrows():
        chrom = str(int(row['chromosome'])) if not pd.isna(row['chromosome']) else None
        if chrom:
            centromere_map[chrom] = {
                'start': row['start'],
                'end': row['end']
            }
    
    colors = {'TAG': '#e74c3c', 'WGD': '#2ecc71'}
    
    # Get chromosomes with data, plot first 5
    chromosomes_to_plot = sorted(subset['chromosome'].unique())[:5]
    
    for chrom in chromosomes_to_plot:
        chr_data = subset[subset['chromosome'] == chrom].copy()
        
        if len(chr_data) == 0:
            continue
        
        print(f"  Chromosome {chrom}: {len(chr_data)} pairs")
        print(f"    Distance range: {chr_data['min_dist_to_end'].min()/1e6:.2f} - {chr_data['min_dist_to_end'].max()/1e6:.2f} Mb")
        
        # Create histogram
        fig, ax = plt.subplots(figsize=(11, 6))
        for dtype in ['TAG', 'WGD']:
            data = chr_data[chr_data['duplication_type'] == dtype]['min_dist_to_end'].dropna()
            if len(data) > 0:
                ax.hist(data / 1_000_000, bins=40, alpha=0.6, label=f'{dtype} (n={len(data)})',
                       color=colors.get(dtype, 'gray'), edgecolor='black')
        
        # Highlight centromere region if available
        if chrom in centromere_map:
            cent = centromere_map[chrom]
            chr_len = chr_lengths.get(chrom, 1)
            
            # Distance from centromere center to nearest end
            cent_center = (cent['start'] + cent['end']) / 2
            dist_to_left = cent_center
            dist_to_right = chr_len - cent_center
            cent_dist_from_end = min(dist_to_left, dist_to_right) / 1_000_000
            
            # Show centromere region
            ax.axvspan(max(0, cent_dist_from_end - PERICENTROMERIC_DISTANCE / 1_000_000),
                      cent_dist_from_end + PERICENTROMERIC_DISTANCE / 1_000_000,
                      alpha=0.2, color='#9b59b6', label='Centromere ±2Mb')
            ax.axvline(cent_dist_from_end, color='#9b59b6', linestyle=':', 
                      linewidth=2, label=f'Centromere ({cent_dist_from_end:.0f}Mb from end)')
        
        # Telomere threshold
        ax.axvline(TELOMERE_DISTANCE / 1_000_000, color='red', linestyle='--',
                  linewidth=2, label=f'Telomere threshold ({TELOMERE_DISTANCE/1_000_000:.0f}Mb)')
        
        ax.set_xlabel('Distance to chromosome end (Mb)', fontsize=12)
        ax.set_ylabel('Frequency', fontsize=12)
        ax.set_title(f'Chromosome {chrom}: Distance to Chromosome Ends\n(Purple = centromere region)', 
                    fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(outdir / f'distance_to_end_chr{chrom}_histogram.png', dpi=300, bbox_inches='tight')
        plt.savefig(outdir / f'distance_to_end_chr{chrom}_histogram.pdf', bbox_inches='tight')
        plt.close()
    
    print(f"  ✓ Saved per-chromosome distance-to-end plots for {len(chromosomes_to_plot)} chromosomes")


def plot_distance_to_end(features_df, outdir: Path):
    """Detailed analysis of distance to chromosome ends - separate plots."""
    print("Generating distance-to-end analysis plots...")
    
    subset = features_df[features_df['same_chromosome']].copy()
    subset = subset[subset['duplication_type'].isin(['TAG', 'WGD'])]
    colors = {'TAG': '#e74c3c', 'WGD': '#2ecc71'}
    
    # Calculate expected centromere distance (roughly in the middle of chromosomes)
    # For Glycine max (~120Mb avg chromosome), centromere is ~40-60Mb from nearest end
    centromere_min_dist = 30_000_000  # 30Mb (conservative estimate from nearest end)
    centromere_max_dist = 60_000_000  # 60Mb (conservative estimate from nearest end)
    
    # 1. Histogram: Distance to end
    fig, ax = plt.subplots(figsize=(10, 6))
    for dtype in ['TAG', 'WGD']:
        data = subset[subset['duplication_type'] == dtype]['min_dist_to_end'].dropna()
        if len(data) > 0:
            ax.hist(data / 1_000_000, bins=60, alpha=0.6, label=dtype,
                   color=colors.get(dtype, 'gray'), edgecolor='black')
    
    # Highlight expected centromere region
    # ax.axvspan(centromere_min_dist / 1_000_000, centromere_max_dist / 1_000_000, 
    #           alpha=0.2, color='#9b59b6', label='Expected centromere region')
    
    # Telomere threshold line
    ax.axvline(TELOMERE_DISTANCE / 1_000_000, color='red', linestyle='--',
              linewidth=2.5, label=f'Telomere threshold ({TELOMERE_DISTANCE/1_000_000:.0f}Mb)')
    
    ax.set_xlabel('Distance to chromosome end (Mb)', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distance to Chromosome Ends: Histogram\n(Purple = expected centromere region)', 
                fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outdir / 'distance_to_end_histogram.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'distance_to_end_histogram.pdf', bbox_inches='tight')
    plt.close()
    
    # 2. CDF
    fig, ax = plt.subplots(figsize=(10, 6))
    for dtype in ['TAG', 'WGD']:
        data = subset[subset['duplication_type'] == dtype]['min_dist_to_end'].dropna() / 1_000_000
        data_sorted = np.sort(data)
        ax.plot(data_sorted, np.arange(1, len(data_sorted) + 1) / len(data_sorted),
               label=dtype, linewidth=2.5, color=colors.get(dtype, 'gray'))
    
    ax.axvline(TELOMERE_DISTANCE / 1_000_000, color='red', linestyle='--',
              linewidth=2, label=f'Threshold ({TELOMERE_DISTANCE/1_000_000:.0f}Mb)')
    ax.set_xlabel('Distance to chromosome end (Mb)', fontsize=12)
    ax.set_ylabel('Cumulative Fraction', fontsize=12)
    ax.set_title('Distance to Chromosome Ends: CDF', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(outdir / 'distance_to_end_cdf.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'distance_to_end_cdf.pdf', bbox_inches='tight')
    plt.close()
    
    # 3. Box plot
    fig, ax = plt.subplots(figsize=(10, 6))
    bp_data = []
    bp_labels = []
    for dtype in ['TAG', 'WGD']:
        data = subset[subset['duplication_type'] == dtype]['min_dist_to_end'].dropna() / 1_000_000
        if len(data) > 0:
            bp_data.append(data.values)
            bp_labels.append(dtype)
    
    bp = ax.boxplot(bp_data, labels=bp_labels, patch_artist=True)
    for patch, dtype in zip(bp['boxes'], bp_labels):
        patch.set_facecolor(colors.get(dtype, 'gray'))
        patch.set_alpha(0.7)
    
    ax.axhline(TELOMERE_DISTANCE / 1_000_000, color='red', linestyle='--',
              linewidth=2, label='Threshold')
    ax.set_ylabel('Distance to chromosome end (Mb)', fontsize=12)
    ax.set_title('Distance to Chromosome Ends: Box Plot', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(outdir / 'distance_to_end_boxplot.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'distance_to_end_boxplot.pdf', bbox_inches='tight')
    plt.close()
    
    # 4. Percentage within threshold
    fig, ax = plt.subplots(figsize=(10, 6))
    pct_within = []
    labels = []
    for dtype in ['TAG', 'WGD']:
        data = subset[subset['duplication_type'] == dtype]
        within = (data['close_to_end']).sum()
        total = len(data)
        pct = within / total * 100 if total > 0 else 0
        pct_within.append(pct)
        labels.append(f'{dtype}\n(n={int(within)}/{total})')
    
    bars = ax.bar(range(len(labels)), pct_within, 
                  color=[colors.get(l.split('\n')[0], 'gray') for l in labels],
                  alpha=0.8, edgecolor='black', width=0.6)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel('Percentage within 4Mb of end (%)', fontsize=12)
    ax.set_title('Pairs Near Chromosome Ends', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    ax.set_ylim(0, 100)
    
    for bar, pct in zip(bars, pct_within):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
               f'{pct:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(outdir / 'distance_to_end_percentage.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'distance_to_end_percentage.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved distance-to-end plots")


def create_summary_report(features_df, outdir: Path):
    """Create comprehensive summary statistics report."""
    print("Generating summary report...")
    
    report_lines = [
        "CHROMOSOMAL STRUCTURAL PATTERNS ANALYSIS",
        "=" * 80,
        ""
    ]
    
    # Overall statistics
    subset = features_df[features_df['duplication_type'].isin(['TAG', 'WGD'])]
    
    report_lines.extend([
        "1. OVERALL SUMMARY",
        f"   Total gene pairs analyzed: {len(subset):,}",
        f"   - TAG duplicates: {(subset['duplication_type'] == 'TAG').sum():,}",
        f"   - WGD duplicates: {(subset['duplication_type'] == 'WGD').sum():,}",
        ""
    ])
    
    # Telomere association
    report_lines.extend([
        "2. TELOMERIC REGION ASSOCIATION",
        "   Gene 1 in telomeric regions:"
    ])
    
    telomere_counts = pd.crosstab(subset['duplication_type'], 
                                  subset['gene1_telomere_region'])
    for dtype in ['TAG', 'WGD']:
        if dtype in telomere_counts.index:
            telomeric = telomere_counts.loc[dtype, 'Telomere-proximal'] if 'Telomere-proximal' in telomere_counts.columns else 0
            total = telomere_counts.loc[dtype].sum()
            pct = telomeric / total * 100 if total > 0 else 0
            report_lines.append(f"   - {dtype}: {int(telomeric)} / {total} ({pct:.1f}%)")
    report_lines.append("")
    
    # Zone distribution
    report_lines.extend([
        "3. CHROMOSOME ZONE DISTRIBUTION (Quintiles)",
        "   Percentages by duplication type:"
    ])
    
    zone_order = ['Q1 (0-20%)', 'Q2 (20-40%)', 'Q3 (40-60%)', 'Q4 (60-80%)', 'Q5 (80-100%)']
    zone_pct = pd.crosstab(subset['duplication_type'], subset['gene1_zone'], normalize='index') * 100
    
    for dtype in ['TAG', 'WGD']:
        if dtype in zone_pct.index:
            report_lines.append(f"\n   {dtype}:")
            for zone in zone_order:
                if zone in zone_pct.columns:
                    pct = zone_pct.loc[dtype, zone]
                    report_lines.append(f"     {zone}: {pct:.1f}%")
    report_lines.append("")
    
    # Ks by location
    report_lines.extend([
        "4. KS VALUES BY CHROMOSOMAL LOCATION",
        "   Mean Ks by region:"
    ])
    
    for dtype in ['TAG', 'WGD']:
        subset_dtype = subset[subset['duplication_type'] == dtype]
        telomeric = subset_dtype[subset_dtype['gene1_telomere_region'] == 'Telomere-proximal']['ks']
        internal = subset_dtype[subset_dtype['gene1_telomere_region'] == 'Internal']['ks']
        
        report_lines.append(f"\n   {dtype}:")
        if len(telomeric) > 0:
            report_lines.append(f"     Telomeric: mean={telomeric.mean():.4f}, median={telomeric.median():.4f}, n={len(telomeric)}")
        if len(internal) > 0:
            report_lines.append(f"     Internal: mean={internal.mean():.4f}, median={internal.median():.4f}, n={len(internal)}")
    report_lines.append("")
    
    # Distance to chromosome ends (same chr only)
    report_lines.extend([
        "5. DISTANCE TO CHROMOSOME ENDS (Same chromosome pairs only)",
        "   Pairs within <4Mb> of chromosome ends:"
    ])
    
    same_chr = subset[subset['same_chromosome']]
    for dtype in ['TAG', 'WGD']:
        subset_dtype = same_chr[same_chr['duplication_type'] == dtype]
        close = subset_dtype['close_to_end'].sum()
        total = len(subset_dtype)
        pct = close / total * 100 if total > 0 else 0
        report_lines.append(f"   - {dtype}: {close} / {total} ({pct:.1f}%)")
    
    report_lines.extend(["", "=" * 80, ""])
    
    # Save report
    report_text = '\n'.join(report_lines)
    report_path = outdir / 'structural_patterns_summary.txt'
    
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    print(f"  ✓ Saved summary report to {report_path}")
    print("\n" + report_text)
    
    return report_text


def main():
    parser = argparse.ArgumentParser(
        description='Analyze correlation between WGD/TAG duplicates and chromosomal structural patterns'
    )
    parser.add_argument('--classified', type=Path, 
                       default=Path('analysis/duplication_types/gene_pairs_classified.tsv'),
                       help='Classified gene pairs TSV')
    parser.add_argument('--protein-info', type=Path, 
                       default=Path('data/protein_info_longest.csv'),
                       help='Protein/gene metadata CSV')
    parser.add_argument('--centromeres', type=Path, default=None,
                       help='(Optional) Centromere locations TSV with columns: chromosome, start, end')
    parser.add_argument('--outdir', type=Path, 
                       default=Path('analysis/structural_analysis'),
                       help='Output directory')
    
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("CHROMOSOMAL STRUCTURAL PATTERN ANALYSIS")
    print("WGD/TAG duplicates vs chromosome architecture (4Mb thresholds)")
    print("=" * 80)
    
    # Load data
    classified_df, protein_df, centromere_df = load_data(
        args.classified, args.protein_info, args.centromeres
    )
    
    # Get chromosome lengths
    chr_lengths = get_chromosome_lengths(protein_df)
    
    # Add structural features
    features_df = add_structural_features(classified_df, protein_df, chr_lengths, centromere_df)
    
    # Save annotated dataset
    features_path = args.outdir / 'gene_pairs_with_structural_features.tsv'
    features_df.to_csv(features_path, sep='\t', index=False, float_format='%.6f')
    print(f"\n✓ Saved annotated dataset to: {features_path}")
    
    # Analyze patterns
    analyze_telomere_association(features_df, args.outdir)
    analyze_zone_distribution(features_df, args.outdir)
    analyze_ks_by_location(features_df, args.outdir)
    
    # Generate visualizations
    plot_telomere_patterns(features_df, args.outdir)
    plot_pericentromeric_patterns(features_df, args.outdir)
    plot_zone_distribution(features_df, args.outdir)
    plot_distance_to_end(features_df, args.outdir)
    plot_distance_to_end_per_chromosome(features_df, classified_df, protein_df, centromere_df, args.outdir)
    
    # Create summary report
    create_summary_report(features_df, args.outdir)
    
    print("\n" + "=" * 80)
    print("✓ ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nAll outputs saved to: {args.outdir}")
    print("\nGenerated files:")
    print("  - gene_pairs_with_structural_features.tsv")
    print("  - structural_telomere_patterns.png/pdf")
    print("  - structural_pericentromeric_patterns.png/pdf (if centromere data provided)")
    print("  - structural_zone_distribution.png/pdf")
    print("  - structural_distance_to_end.png/pdf")
    print("  - structural_patterns_summary.txt")
    print("=" * 80)


if __name__ == '__main__':
    main()
