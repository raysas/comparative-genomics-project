#!/usr/bin/env python3
"""
Plot Ks distribution and related statistics from ks_results_filtered.tsv
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import sys

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

def load_data(tsv_file, max_ks=5.0, sample_size=None):
    """Load and filter Ks results"""
    print(f"Loading {tsv_file}...")
    df = pd.read_csv(tsv_file, sep='\t')
    
    print(f"Total records: {len(df)}")
    print(f"Columns: {list(df.columns)}")
    
    # Filter by Ks value
    df_filtered = df[(df['ks'] >= 0) & (df['ks'] <= max_ks)].copy()
    print(f"After Ks filtering (0 < Ks <= {max_ks}): {len(df_filtered)}")
    
    # Sample if dataset is too large
    if sample_size and len(df_filtered) > sample_size:
        df_filtered = df_filtered.sample(n=sample_size, random_state=42)
        print(f"Sampled to {sample_size} records for faster plotting")
    
    return df_filtered

def plot_ks_distribution(df, output_dir):
    """Create Ks distribution plots"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Ks (Synonymous Substitution Rate) Distribution', fontsize=16, fontweight='bold')
    
    # Histogram
    axes[0, 0].hist(df['ks'], bins=100, edgecolor='black', alpha=0.7, color='steelblue')
    axes[0, 0].set_xlabel('Ks', fontsize=11)
    axes[0, 0].set_ylabel('Frequency', fontsize=11)
    axes[0, 0].set_title('Ks Histogram (100 bins)')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Density plot
    df['ks'].plot(kind='density', ax=axes[0, 1], linewidth=2, color='darkgreen')
    axes[0, 1].set_xlabel('Ks', fontsize=11)
    axes[0, 1].set_title('Ks Density Plot')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Box plot
    axes[1, 0].boxplot(df['ks'], vert=True)
    axes[1, 0].set_ylabel('Ks', fontsize=11)
    axes[1, 0].set_title('Ks Box Plot')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Log-scale histogram
    axes[1, 1].hist(df['ks'], bins=100, edgecolor='black', alpha=0.7, color='coral')
    axes[1, 1].set_yscale('log')
    axes[1, 1].set_xlabel('Ks', fontsize=11)
    axes[1, 1].set_ylabel('Frequency (log scale)', fontsize=11)
    axes[1, 1].set_title('Ks Histogram (Log Scale)')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = Path(output_dir) / 'ks_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_ka_ks_scatter(df, output_dir):
    """Create Ka vs Ks scatter plot"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Scatter plot with density coloring
    scatter = ax.scatter(df['ks'], df['ka'], alpha=0.5, c=df['ka_ks'], 
                        cmap='viridis', s=20, edgecolors='none')
    
    ax.set_xlabel('Ks (Synonymous Substitution Rate)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Ka (Non-synonymous Substitution Rate)', fontsize=12, fontweight='bold')
    ax.set_title('Ka vs Ks Scatter Plot (colored by Ka/Ks ratio)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Ka/Ks Ratio', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    output_file = Path(output_dir) / 'ka_vs_ks_scatter.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_ka_ks_distribution(df, output_dir):
    """Create Ka/Ks distribution plots"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Ka/Ks Ratio Distribution', fontsize=16, fontweight='bold')
    
    # Histogram
    axes[0].hist(df['ka_ks'], bins=100, edgecolor='black', alpha=0.7, color='mediumseagreen')
    axes[0].set_xlabel('Ka/Ks', fontsize=11)
    axes[0].set_ylabel('Frequency', fontsize=11)
    axes[0].set_title('Ka/Ks Histogram')
    axes[0].axvline(1.0, color='red', linestyle='--', linewidth=2, label='Ka/Ks = 1')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Density plot
    df['ka_ks'].plot(kind='density', ax=axes[1], linewidth=2, color='darkred')
    axes[1].axvline(1.0, color='red', linestyle='--', linewidth=2, label='Ka/Ks = 1')
    axes[1].set_xlabel('Ka/Ks', fontsize=11)
    axes[1].set_title('Ka/Ks Density Plot')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = Path(output_dir) / 'ka_ks_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_family_distribution(df, output_dir, top_n=20):
    """Create family-wise Ks distribution"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('Family-wise Ks Distribution', fontsize=16, fontweight='bold')
    
    # Top families by count
    top_families = df['family'].value_counts().head(top_n)
    family_ks_data = [df[df['family'] == fam]['ks'].values for fam in top_families.index]
    
    # Box plot
    bp = axes[0].boxplot(family_ks_data, labels=top_families.index, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    axes[0].set_xlabel('Family', fontsize=11)
    axes[0].set_ylabel('Ks', fontsize=11)
    axes[0].set_title(f'Ks Distribution for Top {top_n} Families (by count)')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # Mean Ks per family
    family_mean_ks = df.groupby('family')['ks'].mean().sort_values(ascending=False).head(top_n)
    axes[1].barh(range(len(family_mean_ks)), family_mean_ks.values, color='steelblue')
    axes[1].set_yticks(range(len(family_mean_ks)))
    axes[1].set_yticklabels(family_mean_ks.index)
    axes[1].set_xlabel('Mean Ks', fontsize=11)
    axes[1].set_title(f'Mean Ks for Top {top_n} Families (by mean Ks)')
    axes[1].grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    output_file = Path(output_dir) / 'family_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def plot_length_vs_ks(df, output_dir):
    """Create length vs Ks scatter plot"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    scatter = ax.scatter(df['length'], df['ks'], alpha=0.4, c=df['ka_ks'], 
                        cmap='plasma', s=20, edgecolors='none')
    
    ax.set_xlabel('Alignment Length (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Ks', fontsize=12, fontweight='bold')
    ax.set_title('Alignment Length vs Ks (colored by Ka/Ks)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Ka/Ks Ratio', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    output_file = Path(output_dir) / 'length_vs_ks.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()

def print_statistics(df):
    """Print summary statistics"""
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"\nKs (Synonymous Substitution Rate):")
    print(f"  Mean:     {df['ks'].mean():.4f}")
    print(f"  Median:   {df['ks'].median():.4f}")
    print(f"  Std Dev:  {df['ks'].std():.4f}")
    print(f"  Min:      {df['ks'].min():.4f}")
    print(f"  Max:      {df['ks'].max():.4f}")
    print(f"  Q1:       {df['ks'].quantile(0.25):.4f}")
    print(f"  Q3:       {df['ks'].quantile(0.75):.4f}")
    
    print(f"\nKa (Non-synonymous Substitution Rate):")
    print(f"  Mean:     {df['ka'].mean():.4f}")
    print(f"  Median:   {df['ka'].median():.4f}")
    print(f"  Std Dev:  {df['ka'].std():.4f}")
    print(f"  Min:      {df['ka'].min():.4f}")
    print(f"  Max:      {df['ka'].max():.4f}")
    
    print(f"\nKa/Ks Ratio:")
    print(f"  Mean:     {df['ka_ks'].mean():.4f}")
    print(f"  Median:   {df['ka_ks'].median():.4f}")
    print(f"  Std Dev:  {df['ka_ks'].std():.4f}")
    print(f"  Min:      {df['ka_ks'].min():.4f}")
    print(f"  Max:      {df['ka_ks'].max():.4f}")
    print(f"  Count (Ka/Ks < 1): {(df['ka_ks'] < 1).sum()} (purifying selection)")
    print(f"  Count (Ka/Ks = 1): {(df['ka_ks'] == 1).sum()} (neutral)")
    print(f"  Count (Ka/Ks > 1): {(df['ka_ks'] > 1).sum()} (positive selection)")
    
    print(f"\nAlignment Length:")
    print(f"  Mean:     {df['length'].mean():.1f} bp")
    print(f"  Median:   {df['length'].median():.1f} bp")
    print(f"  Min:      {df['length'].min():.1f} bp")
    print(f"  Max:      {df['length'].max():.1f} bp")
    
    print(f"\nFamilies: {df['family'].nunique()}")
    print(f"Gene pairs: {len(df)}")
    print("="*60 + "\n")

def main():
    # File path
    tsv_file = Path(__file__).parent / 'ks_results_filtered.tsv'
    output_dir = Path(__file__).parent / 'ks_plots'
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    # Load data (sample to 50k records for speed if needed)
    df = load_data(str(tsv_file), max_ks=5.0, sample_size=50000)
    
    # Print statistics
    print_statistics(df)
    
    # Create plots
    print("Creating plots...")
    plot_ks_distribution(df, output_dir)
    plot_ka_ks_scatter(df, output_dir)
    plot_ka_ks_distribution(df, output_dir)
    plot_family_distribution(df, output_dir, top_n=20)
    plot_length_vs_ks(df, output_dir)
    
    print(f"\n✓ All plots saved to: {output_dir}")

if __name__ == '__main__':
    main()
