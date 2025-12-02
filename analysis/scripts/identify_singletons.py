#!/usr/bin/env python3
"""
Identify singleton genes (non-duplicated genes) and analyze their characteristics.

Singletons are genes that:
- Are not part of any gene family
- Have no duplicates based on BLAST similarity

Usage:
    python3 identify_singletons.py

Outputs:
    - singletons.tsv: List of singleton genes with metadata
    - singletons_summary.tsv: Summary statistics
    - singletons_plots.png: Visualizations
"""
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10

# File paths
PROTEIN_INFO = "data/protein_info_longest.csv"
FAMILY_FILE = "output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv"
if len(sys.argv) > 1:
    FAMILY_FILE = sys.argv[1]

OUTPUT_DIR = "output/singletons/"

# Create output directory
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

def load_data():
    """Load protein info and family assignments."""
    print("Loading data...")
    
    # Load all genes
    protein_df = pd.read_csv(PROTEIN_INFO)
    print(f"  ✓ Loaded {len(protein_df)} total genes")
    
    # Load gene families
    families_df = pd.read_csv(FAMILY_FILE, sep='\t')
    print(f"  ✓ Loaded {len(families_df)} genes in families")
    
    return protein_df, families_df

def identify_singletons(protein_df, families_df):
    """Identify singleton genes."""
    print("\nIdentifying singletons...")
    
    # Get duplicated gene IDs
    duplicated_genes = set(families_df['geneName'])
    
    # Mark singletons
    protein_df['is_singleton'] = ~protein_df['peptide_id'].isin(duplicated_genes)
    
    # Extract singletons
    singletons_df = protein_df[protein_df['is_singleton']].copy()
    
    print(f"  Total genes: {len(protein_df)}")
    print(f"  Duplicated genes: {len(duplicated_genes)}")
    print(f"  Singleton genes: {len(singletons_df)}")
    print(f"  Singleton percentage: {len(singletons_df)/len(protein_df)*100:.2f}%")
    
    return singletons_df, protein_df

def analyze_by_chromosome(singletons_df, protein_df):
    """Analyze singleton distribution by chromosome."""
    print("\nAnalyzing chromosome distribution...")
    
    # Count by chromosome
    chr_all = protein_df['chromosome'].value_counts().sort_index()
    chr_singletons = singletons_df['chromosome'].value_counts().sort_index()
    chr_duplicated = chr_all - chr_singletons
    
    chr_stats = pd.DataFrame({
        'Total_genes': chr_all,
        'Singletons': chr_singletons,
        'Duplicated': chr_duplicated,
        'Pct_Singleton': (chr_singletons / chr_all * 100).round(2)
    }).fillna(0)
    
    # Save
    chr_stats.to_csv(f"{OUTPUT_DIR}/singletons_by_chromosome.tsv", sep='\t')
    print(f"  ✓ Saved chromosome statistics")
    
    return chr_stats

def create_summary_statistics(singletons_df, protein_df):
    """Create summary statistics."""
    print("\nGenerating summary statistics...")
    
    total_genes = len(protein_df)
    n_singletons = len(singletons_df)
    n_duplicated = total_genes - n_singletons
    
    # Length statistics
    singleton_lengths = singletons_df['length'].describe()
    duplicated_lengths = protein_df[~protein_df['is_singleton']]['length'].describe()
    
    summary = {
        'Metric': [
            'Total genes',
            'Singleton genes',
            'Duplicated genes',
            'Singleton percentage',
            'Mean length (singletons)',
            'Mean length (duplicated)',
            'Median length (singletons)',
            'Median length (duplicated)'
        ],
        'Value': [
            total_genes,
            n_singletons,
            n_duplicated,
            f"{n_singletons/total_genes*100:.2f}%",
            f"{singleton_lengths['mean']:.1f}",
            f"{duplicated_lengths['mean']:.1f}",
            f"{singleton_lengths['50%']:.1f}",
            f"{duplicated_lengths['50%']:.1f}"
        ]
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(f"{OUTPUT_DIR}/singletons_summary.tsv", sep='\t', index=False)
    print(f"  ✓ Saved summary statistics")
    
    return summary_df

def plot_visualizations(singletons_df, protein_df, chr_stats):
    """Create comprehensive visualizations."""
    print("\nGenerating visualizations...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    fig.suptitle('Singleton Gene Analysis', fontsize=16, fontweight='bold')
    
    # 1. Pie chart: Singleton vs Duplicated
    ax1 = fig.add_subplot(gs[0, 0])
    counts = [len(singletons_df), len(protein_df) - len(singletons_df)]
    labels = ['Singletons', 'Duplicated']
    colors = ['lightcoral', 'steelblue']
    wedges, texts, autotexts = ax1.pie(counts, labels=labels, autopct='%1.1f%%',
                                         colors=colors, startangle=90)
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    ax1.set_title('Genome Composition', fontweight='bold')
    
    # 2. Bar chart: Count comparison
    ax2 = fig.add_subplot(gs[0, 1:])
    x = ['Singletons', 'Duplicated']
    y = counts
    bars = ax2.bar(x, y, color=colors, edgecolor='black', linewidth=1.5)
    ax2.set_ylabel('Number of genes')
    ax2.set_title('Gene Count Comparison', fontweight='bold')
    ax2.grid(alpha=0.3, axis='y')
    for bar, val in zip(bars, y):
        ax2.text(bar.get_x() + bar.get_width()/2, val + max(y)*0.02,
                f'{val:,}', ha='center', va='bottom', fontweight='bold', fontsize=11)
    
    # 3. Chromosome distribution (stacked bar)
    ax3 = fig.add_subplot(gs[1, :])
    chr_stats[['Singletons', 'Duplicated']].plot(kind='bar', stacked=True, ax=ax3,
                                                   color=['lightcoral', 'steelblue'])
    ax3.set_xlabel('Chromosome')
    ax3.set_ylabel('Number of genes')
    ax3.set_title('Gene Distribution by Chromosome', fontweight='bold')
    ax3.legend(['Singletons', 'Duplicated'])
    ax3.grid(alpha=0.3, axis='y')
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=0)
    
    # 4. Percentage of singletons per chromosome
    ax4 = fig.add_subplot(gs[2, 0])
    chr_stats['Pct_Singleton'].plot(kind='bar', ax=ax4, color='coral')
    ax4.set_xlabel('Chromosome')
    ax4.set_ylabel('% Singletons')
    ax4.set_title('Singleton Percentage by Chromosome', fontweight='bold')
    ax4.axhline(chr_stats['Pct_Singleton'].mean(), color='red', linestyle='--',
                label=f"Mean: {chr_stats['Pct_Singleton'].mean():.1f}%")
    ax4.legend()
    ax4.grid(alpha=0.3, axis='y')
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # 5. Length distribution comparison
    ax5 = fig.add_subplot(gs[2, 1])
    singleton_lengths = singletons_df['length']
    duplicated_lengths = protein_df[~protein_df['is_singleton']]['length']
    
    ax5.hist([singleton_lengths, duplicated_lengths], bins=50, label=['Singletons', 'Duplicated'],
            color=['lightcoral', 'steelblue'], alpha=0.7, edgecolor='black', linewidth=0.5)
    ax5.set_xlabel('Protein length (aa)')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Protein Length Distribution', fontweight='bold')
    ax5.legend()
    ax5.grid(alpha=0.3, axis='y')
    
    # 6. Box plot: Length comparison
    ax6 = fig.add_subplot(gs[2, 2])
    data_to_plot = [singleton_lengths, duplicated_lengths]
    bp = ax6.boxplot(data_to_plot, labels=['Singletons', 'Duplicated'], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightcoral')
    bp['boxes'][1].set_facecolor('steelblue')
    ax6.set_ylabel('Protein length (aa)')
    ax6.set_title('Length Comparison', fontweight='bold')
    ax6.grid(alpha=0.3, axis='y')
    
    # Add mean markers
    means = [singleton_lengths.mean(), duplicated_lengths.mean()]
    ax6.plot([1, 2], means, 'r*', markersize=12, label='Mean')
    ax6.legend()
    
    plt.savefig(f"{OUTPUT_DIR}/singletons_analysis.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{OUTPUT_DIR}/singletons_analysis.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved visualizations")

def main():
    print("="*60)
    print("SINGLETON GENE IDENTIFICATION")
    print("="*60)
    
    # Load data
    protein_df, families_df = load_data()
    
    # Identify singletons
    singletons_df, protein_df = identify_singletons(protein_df, families_df)
    
    # Save singleton list
    output_file = f"{OUTPUT_DIR}/singletons.tsv"
    singletons_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n  ✓ Saved singleton list to: {output_file}")
    
    # Analyze by chromosome
    chr_stats = analyze_by_chromosome(singletons_df, protein_df)
    
    # Summary statistics
    summary_df = create_summary_statistics(singletons_df, protein_df)
    
    # Visualizations
    plot_visualizations(singletons_df, protein_df, chr_stats)
    
    print("\n" + "="*60)
    print("✓ SINGLETON ANALYSIS COMPLETE")
    print("="*60)
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
