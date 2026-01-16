#!/usr/bin/env python3
"""
Plot orientation of TAG gene pairs from TAGs_genecount_10genes.tsv

Creates:
 - bar plot of Same vs Opposite orientation
 - pie chart of same vs opposite orientation
 - distance and gene separation plots by orientation
 - saves summary counts TSV

Usage:
    python3 plot_TAG_orientation.py \
        --input analysis/duplication_types/TAGs/TAGs_genecount_10genes.tsv \
        --outdir analysis/duplication_types/TAGs
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency

sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 150


def load_tags(input_path: Path):
    df = pd.read_csv(input_path, sep='\t')
    return df


def categorize_orientation(df: pd.DataFrame) -> pd.DataFrame:
    """Add orientation category as Same or Opposite."""
    # Use Same_Orientation column if available
    if 'Same_Orientation' in df.columns:
        df['Orientation'] = df['Same_Orientation'].astype(str).map(
            lambda x: 'Same' if x.lower() in ['true', '1', 't', 'yes'] else 'Opposite'
        )
    else:
        # Otherwise derive from strands
        df['Strand1'] = df['Strand1'].astype(str).str.strip()
        df['Strand2'] = df['Strand2'].astype(str).str.strip()
        df['Orientation'] = (df['Strand1'] == df['Strand2']).map(
            lambda x: 'Same' if x else 'Opposite'
        )
    
    return df


def plot_orientation_bar(df: pd.DataFrame, outdir: Path):
    """Bar plot showing Same vs Opposite orientations with chi-square test."""
    outdir.mkdir(parents=True, exist_ok=True)

    # Count orientations
    ori_counts = df['Orientation'].value_counts()
    
    # Ensure both categories present in order
    ori_order = ['Same', 'Opposite']
    ori_counts = ori_counts.reindex(ori_order, fill_value=0)

    # Chi-square test (null: equal proportions)
    chi2, p_val, dof, expected = chi2_contingency([ori_counts.values])

    # Bar plot
    fig, ax = plt.subplots(figsize=(9, 7))
    colors = {'Same': '#2ecc71', 'Opposite': '#e74c3c'}
    color_list = [colors[ori] for ori in ori_counts.index]
    
    bars = ax.bar(ori_counts.index, ori_counts.values, color=color_list, 
                  edgecolor='black', alpha=0.85, linewidth=2)
    
    ax.set_ylabel('Number of TAG pairs', fontsize=12, fontweight='bold')
    ax.set_xlabel('Orientation', fontsize=12, fontweight='bold')
    ax.set_title('TAG Orientation Distribution', fontsize=14, fontweight='bold', pad=25)
    ax.grid(alpha=0.3, axis='y')

    # Add count and percentage labels
    total = len(df)
    for bar, ori in zip(bars, ori_counts.index):
        count = ori_counts[ori]
        pct = count / total * 100
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, height, 
                f'{int(count):,}\n({pct:.1f}%)', 
                ha='center', va='bottom', fontweight='bold', fontsize=11)

    # Add chi-square test result as text
    sig_marker = '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
    test_text = f'χ² = {chi2:.4f}, p = {p_val:.4e} {sig_marker}'
    ax.text(0.5, -0.15, test_text, transform=ax.transAxes, 
            ha='center', fontsize=10, style='italic',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.tight_layout()
    out_png = outdir / 'TAG_orientation_bar.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(outdir / 'TAG_orientation_bar.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ Saved bar plot to {out_png}")
    print(f"    Chi-square test: {test_text}")


def plot_orientation_pie(df: pd.DataFrame, outdir: Path):
    """Pie chart showing Same vs Opposite orientations."""
    ori_counts = df['Orientation'].value_counts()
    ori_counts = ori_counts.reindex(['Same', 'Opposite'], fill_value=0)

    fig, ax = plt.subplots(figsize=(8, 7))
    colors = ['#2ecc71', '#e74c3c']
    
    wedges, texts, autotexts = ax.pie(
        ori_counts.values, 
        labels=[f'{ori}\n(n={int(ori_counts[ori]):,})' for ori in ori_counts.index],
        autopct='%1.1f%%', 
        colors=colors, 
        startangle=90,
        textprops={'fontsize': 11, 'fontweight': 'bold'},
        explode=(0.05, 0.05)
    )
    
    ax.set_title('TAG Orientation Distribution', fontsize=14, fontweight='bold', pad=25)
    
    fig.tight_layout()
    out_png = outdir / 'TAG_orientation_pie.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(outdir / 'TAG_orientation_pie.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ Saved pie chart to {out_png}")


def plot_orientation_by_distance(df: pd.DataFrame, outdir: Path):
    """Plots comparing distance/separation by orientation."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('TAG Orientation vs Gene Separation/Distance', fontsize=14, fontweight='bold', y=0.995)
    
    colors = {'Same': '#2ecc71', 'Opposite': '#e74c3c'}
    
    # 1. Histogram: Genes Apart
    ax = axes[0, 0]
    if 'Genes_Apart' in df.columns:
        for ori in ['Same', 'Opposite']:
            data = df[df['Orientation'] == ori]['Genes_Apart']
            if len(data) > 0:
                ax.hist(data, bins=range(1, min(int(data.max())+2, 51)), 
                       alpha=0.6, label=ori, color=colors[ori], edgecolor='black')
        ax.set_xlabel('Genes Apart', fontsize=11, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax.set_title('Gene Separation Distribution', pad=15)
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3, axis='y')
        ax.set_xlim(0, 50)
    
    # 2. Box plot: Genes Apart
    ax = axes[0, 1]
    if 'Genes_Apart' in df.columns:
        bp_data = [df[df['Orientation'] == ori]['Genes_Apart'].values 
                   for ori in ['Same', 'Opposite']]
        bp = ax.boxplot(bp_data, labels=['Same', 'Opposite'], patch_artist=True)
        for patch, ori in zip(bp['boxes'], ['Same', 'Opposite']):
            patch.set_facecolor(colors[ori])
            patch.set_alpha(0.7)
        ax.set_ylabel('Genes Apart', fontsize=11, fontweight='bold')
        ax.set_title('Gene Separation (Box Plot)', pad=15)
        ax.grid(alpha=0.3, axis='y')
    
    # 3. Histogram: Distance in Mb
    ax = axes[1, 0]
    if 'Distance_bp' in df.columns:
        for ori in ['Same', 'Opposite']:
            data = df[df['Orientation'] == ori]['Distance_bp'] / 1_000_000
            if len(data) > 0:
                ax.hist(data, bins=50, alpha=0.6, label=ori, 
                       color=colors[ori], edgecolor='black')
        ax.set_xlabel('Distance (Mb)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=11, fontweight='bold')
        ax.set_title('Base Pair Distance Distribution', pad=15)
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3, axis='y')
    
    # 4. Box plot: Distance
    ax = axes[1, 1]
    if 'Distance_bp' in df.columns:
        bp_data = [df[df['Orientation'] == ori]['Distance_bp'].values / 1_000_000
                   for ori in ['Same', 'Opposite']]
        bp = ax.boxplot(bp_data, labels=['Same', 'Opposite'], patch_artist=True)
        for patch, ori in zip(bp['boxes'], ['Same', 'Opposite']):
            patch.set_facecolor(colors[ori])
            patch.set_alpha(0.7)
        ax.set_ylabel('Distance (Mb)', fontsize=11, fontweight='bold')
        ax.set_title('Base Pair Distance (Box Plot)', pad=15)
        ax.grid(alpha=0.3, axis='y')
    
    fig.tight_layout()
    out_png = outdir / 'TAG_orientation_vs_distance.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(outdir / 'TAG_orientation_vs_distance.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ Saved distance plots to {out_png}")


def save_summary(df: pd.DataFrame, outdir: Path):
    """Save summary statistics."""
    summary_rows = []
    
    for ori in ['Same', 'Opposite']:
        subset = df[df['Orientation'] == ori]
        row = {
            'Orientation': ori,
            'N_pairs': len(subset),
            'Percentage': f"{len(subset)/len(df)*100:.2f}%",
        }
        
        if 'Genes_Apart' in df.columns:
            row['Mean_Genes_Apart'] = subset['Genes_Apart'].mean()
            row['Median_Genes_Apart'] = subset['Genes_Apart'].median()
        
        if 'Distance_bp' in df.columns:
            row['Mean_Distance_Mb'] = subset['Distance_bp'].mean() / 1_000_000
            row['Median_Distance_Mb'] = subset['Distance_bp'].median() / 1_000_000
        
        summary_rows.append(row)
    
    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / 'TAG_orientation_summary.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False, float_format='%.4f')
    
    print(f"\n  ✓ Saved summary to {summary_path}")
    print("\nSummary Statistics:")
    print(summary_df.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description='Plot TAG orientation patterns (Same vs Opposite)')
    parser.add_argument('--input', type=Path, required=True, help='TAGs TSV input file')
    parser.add_argument('--outdir', type=Path, default=Path('analysis/duplication_types/TAGs'), 
                       help='Output directory')
    args = parser.parse_args()

    print("=" * 80)
    print("TAG ORIENTATION ANALYSIS")
    print("=" * 80)
    
    df = load_tags(args.input)
    print(f"\nLoaded {len(df):,} TAG pairs")
    
    df = categorize_orientation(df)
    
    # Generate plots
    print("\nGenerating plots...")
    plot_orientation_bar(df, args.outdir)
    plot_orientation_pie(df, args.outdir)
    plot_orientation_by_distance(df, args.outdir)
    
    # Save summary
    save_summary(df, args.outdir)
    
    print("\n" + "=" * 80)
    print("✓ ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nAll outputs saved to: {args.outdir}")
    print("\nGenerated files:")
    print("  - TAG_orientation_bar.png/pdf")
    print("  - TAG_orientation_pie.png/pdf")
    print("  - TAG_orientation_vs_distance.png/pdf")
    print("  - TAG_orientation_summary.tsv")
    print("=" * 80)

if __name__ == '__main__':
    main()
