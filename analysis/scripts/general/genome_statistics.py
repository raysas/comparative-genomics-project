#!/usr/bin/env python3
"""
Comprehensive genome and pipeline statistics analysis.
Generates summary statistics and visualizations for Phase 1-4 results.

Usage:
    python3 genome_statistics.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10

# File paths
PROTEIN_INFO = "data/protein_info_longest.csv"
FAMILY_FILE = "output/filtered/protein_families_filtered_blast_results_id50_qcov90_scov90_wcol12_network.tsv"
BLAST_FILE = "output/blast_output/blast_results_with_coverage.tsv"
KS_FILE = "output/ks_results/ks_results_filtered.tsv"
OUTPUT_DIR = "output/statistics/"

# Create output directory
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

def load_data():
    """Load all necessary data files."""
    print("Loading data files...")
    
    # Protein info (all genes)
    protein_df = pd.read_csv(PROTEIN_INFO)
    print(f"  ✓ Loaded {len(protein_df)} total genes")
    
    # Gene families
    families_df = pd.read_csv(FAMILY_FILE, sep='\t')
    print(f"  ✓ Loaded {len(families_df)} genes in families")
    
    # BLAST results (if available)
    try:
        blast_df = pd.read_csv(BLAST_FILE, sep='\t')
        print(f"  ✓ Loaded {len(blast_df)} BLAST pairs")
    except:
        blast_df = None
        print(f"  ⚠ BLAST results not found")
    
    # Ks results (if available)
    try:
        ks_df = pd.read_csv(KS_FILE, sep='\t')
        print(f"  ✓ Loaded {len(ks_df)} Ks pairs")
    except:
        ks_df = None
        print(f"  ⚠ Ks results not found")
    
    return protein_df, families_df, blast_df, ks_df

def compute_basic_statistics(protein_df, families_df):
    """Compute basic genome and duplication statistics."""
    print("\n=== Basic Statistics ===")
    
    total_genes = len(protein_df)
    duplicated_genes = len(families_df)
    singleton_genes = total_genes - duplicated_genes
    n_families = families_df['family'].nunique()
    
    # Family sizes
    family_sizes = families_df.groupby('family').size()
    
    stats = {
        'Total genes in genome': total_genes,
        'Duplicated genes': duplicated_genes,
        'Singleton genes': singleton_genes,
        'Percentage duplicated': f"{duplicated_genes/total_genes*100:.2f}%",
        'Number of families': n_families,
        'Mean family size': f"{family_sizes.mean():.2f}",
        'Median family size': int(family_sizes.median()),
        'Max family size': int(family_sizes.max()),
        'Min family size': int(family_sizes.min())
    }
    
    for key, val in stats.items():
        print(f"  {key}: {val}")
    
    # Save to file
    stats_df = pd.DataFrame(list(stats.items()), columns=['Metric', 'Value'])
    stats_df.to_csv(f"{OUTPUT_DIR}/basic_statistics.tsv", sep='\t', index=False)
    
    return stats, family_sizes

def plot_chromosome_distribution(protein_df, families_df):
    """Plot gene distribution across chromosomes."""
    print("\nGenerating chromosome distribution plots...")
    
    # Count genes per chromosome
    chr_counts = protein_df['chromosome'].value_counts().sort_index()
    
    # Count duplicated genes per chromosome
    duplicated_genes = set(families_df['geneName'])
    protein_df['is_duplicated'] = protein_df['peptide_id'].isin(duplicated_genes)
    chr_dup_counts = protein_df[protein_df['is_duplicated']].groupby('chromosome').size()
    
    # Merge counts
    chr_stats = pd.DataFrame({
        'Total': chr_counts,
        'Duplicated': chr_dup_counts
    }).fillna(0)
    chr_stats['Singleton'] = chr_stats['Total'] - chr_stats['Duplicated']
    chr_stats['Pct_Duplicated'] = (chr_stats['Duplicated'] / chr_stats['Total'] * 100).round(2)
    
    # Save stats
    chr_stats.to_csv(f"{OUTPUT_DIR}/genes_per_chromosome.tsv", sep='\t')
    
    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Stacked bar plot
    ax1 = axes[0]
    chr_stats[['Duplicated', 'Singleton']].plot(kind='bar', stacked=True, 
                                                  ax=ax1, color=['steelblue', 'lightgray'])
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Number of genes')
    ax1.set_title('Gene Distribution by Chromosome', fontweight='bold')
    ax1.legend(['Duplicated', 'Singleton'])
    ax1.grid(alpha=0.3, axis='y')
    
    # Percentage plot
    ax2 = axes[1]
    chr_stats['Pct_Duplicated'].plot(kind='bar', ax=ax2, color='coral')
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('% Duplicated')
    ax2.set_title('Percentage of Duplicated Genes by Chromosome', fontweight='bold')
    ax2.grid(alpha=0.3, axis='y')
    ax2.axhline(chr_stats['Pct_Duplicated'].mean(), color='red', linestyle='--', 
                label=f"Mean: {chr_stats['Pct_Duplicated'].mean():.1f}%")
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/genes_per_chromosome.png", dpi=300)
    plt.savefig(f"{OUTPUT_DIR}/genes_per_chromosome.pdf")
    plt.close()
    
    print(f"  ✓ Saved chromosome distribution plots")
    
    return chr_stats

def plot_family_size_distribution(family_sizes):
    """Plot family size distribution."""
    print("\nGenerating family size distribution plots...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Histogram
    ax1 = axes[0]
    ax1.hist(family_sizes, bins=50, edgecolor='black', color='steelblue', alpha=0.7)
    ax1.set_xlabel('Family size (number of genes)')
    ax1.set_ylabel('Number of families')
    ax1.set_title('Distribution of Gene Family Sizes', fontweight='bold')
    ax1.axvline(family_sizes.mean(), color='red', linestyle='--', 
                label=f'Mean: {family_sizes.mean():.1f}')
    ax1.axvline(family_sizes.median(), color='green', linestyle='--', 
                label=f'Median: {family_sizes.median():.0f}')
    ax1.legend()
    ax1.grid(alpha=0.3, axis='y')
    
    # Log scale
    ax2 = axes[1]
    size_counts = family_sizes.value_counts().sort_index()
    ax2.bar(size_counts.index, size_counts.values, edgecolor='black', color='darkgreen', alpha=0.7)
    ax2.set_xlabel('Family size (number of genes)')
    ax2.set_ylabel('Number of families')
    ax2.set_title('Family Size Distribution (counts)', fontweight='bold')
    ax2.set_yscale('log')
    ax2.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/family_size_distribution.png", dpi=300)
    plt.savefig(f"{OUTPUT_DIR}/family_size_distribution.pdf")
    plt.close()
    
    # Save top families
    top_families = family_sizes.nlargest(20)
    top_families_df = pd.DataFrame({
        'Family': top_families.index,
        'Size': top_families.values
    })
    top_families_df.to_csv(f"{OUTPUT_DIR}/top_20_families.tsv", sep='\t', index=False)
    
    print(f"  ✓ Saved family size distribution plots")
    print(f"  ✓ Top 10 largest families:")
    for fam, size in top_families.head(10).items():
        print(f"      Family {fam}: {size} genes")

def plot_pipeline_filtering_effects(blast_df, families_df, ks_df):
    """Show how data flows through the pipeline."""
    print("\nGenerating pipeline filtering visualization...")
    
    # Calculate numbers at each stage
    total_blast_pairs = len(blast_df) if blast_df is not None else 0
    total_family_genes = len(families_df)
    total_family_pairs = len(families_df) * (len(families_df) - 1) / 2  # Approximate
    total_ks_pairs = len(ks_df) if ks_df is not None else 0
    
    # Create funnel plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    stages = ['BLAST pairs\n(all-vs-all)', 'After filtering\n(coverage, identity)', 
              'Gene pairs\n(in families)', 'Ks calculated\n(with alignments)']
    values = [total_blast_pairs, total_blast_pairs, total_family_genes, total_ks_pairs]
    
    # Only include stages with data
    if total_blast_pairs == 0:
        stages = stages[2:]
        values = values[2:]
    if total_ks_pairs == 0:
        stages = stages[:-1]
        values = values[:-1]
    
    colors = plt.cm.Blues(np.linspace(0.4, 0.9, len(stages)))
    y_positions = np.arange(len(stages))
    
    bars = ax.barh(y_positions, values, color=colors, edgecolor='black', linewidth=1.5)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, values)):
        ax.text(val + max(values)*0.02, i, f'{int(val):,}', 
                va='center', fontweight='bold', fontsize=11)
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels(stages)
    ax.set_xlabel('Number of pairs / genes', fontsize=12)
    ax.set_title('Pipeline Filtering Effects', fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/pipeline_filtering.png", dpi=300)
    plt.savefig(f"{OUTPUT_DIR}/pipeline_filtering.pdf")
    plt.close()
    
    print(f"  ✓ Saved pipeline filtering visualization")

def create_summary_figure(stats, chr_stats, family_sizes):
    """Create a comprehensive summary figure."""
    print("\nCreating comprehensive summary figure...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Title
    fig.suptitle('Glycine max Genome Duplication Analysis - Summary Statistics', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # 1. Basic stats text
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.axis('off')
    stats_text = '\n'.join([f"{k}: {v}" for k, v in stats.items()])
    ax1.text(0.1, 0.5, stats_text, fontsize=10, verticalalignment='center',
             family='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax1.set_title('Genome Statistics', fontweight='bold', loc='left')
    
    # 2. Duplication pie chart
    ax2 = fig.add_subplot(gs[0, 1])
    dup_counts = [stats['Duplicated genes'], stats['Singleton genes']]
    colors_pie = ['steelblue', 'lightgray']
    wedges, texts, autotexts = ax2.pie(dup_counts, labels=['Duplicated', 'Singleton'], 
                                         autopct='%1.1f%%', colors=colors_pie, startangle=90)
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    ax2.set_title('Gene Duplication Status', fontweight='bold')
    
    # 3. Family size distribution
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.hist(family_sizes, bins=30, edgecolor='black', color='darkgreen', alpha=0.7)
    ax3.set_xlabel('Family size')
    ax3.set_ylabel('Count')
    ax3.set_title('Family Size Distribution', fontweight='bold')
    ax3.grid(alpha=0.3, axis='y')
    
    # 4. Genes per chromosome
    ax4 = fig.add_subplot(gs[1, :])
    chr_stats[['Duplicated', 'Singleton']].plot(kind='bar', stacked=True, ax=ax4,
                                                  color=['steelblue', 'lightgray'])
    ax4.set_xlabel('Chromosome')
    ax4.set_ylabel('Number of genes')
    ax4.set_title('Gene Distribution Across Chromosomes', fontweight='bold')
    ax4.legend(['Duplicated', 'Singleton'])
    ax4.grid(alpha=0.3, axis='y')
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=0)
    
    # 5. Duplication percentage per chromosome
    ax5 = fig.add_subplot(gs[2, :])
    chr_stats['Pct_Duplicated'].plot(kind='bar', ax=ax5, color='coral')
    ax5.set_xlabel('Chromosome')
    ax5.set_ylabel('% Duplicated')
    ax5.set_title('Percentage of Duplicated Genes by Chromosome', fontweight='bold')
    ax5.axhline(chr_stats['Pct_Duplicated'].mean(), color='red', linestyle='--',
                label=f"Genome mean: {chr_stats['Pct_Duplicated'].mean():.1f}%")
    ax5.legend()
    ax5.grid(alpha=0.3, axis='y')
    plt.setp(ax5.xaxis.get_majorticklabels(), rotation=0)
    
    plt.savefig(f"{OUTPUT_DIR}/summary_figure.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{OUTPUT_DIR}/summary_figure.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved comprehensive summary figure")

def main():
    print("="*60)
    print("GENOME STATISTICS ANALYSIS")
    print("="*60)
    
    # Load data
    protein_df, families_df, blast_df, ks_df = load_data()
    
    # Basic statistics
    stats, family_sizes = compute_basic_statistics(protein_df, families_df)
    
    # Chromosome distribution
    chr_stats = plot_chromosome_distribution(protein_df, families_df)
    
    # Family size distribution
    plot_family_size_distribution(family_sizes)
    
    # Pipeline filtering effects
    plot_pipeline_filtering_effects(blast_df, families_df, ks_df)
    
    # Summary figure
    create_summary_figure(stats, chr_stats, family_sizes)
    
    print("\n" + "="*60)
    print("✓ ANALYSIS COMPLETE")
    print("="*60)
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
