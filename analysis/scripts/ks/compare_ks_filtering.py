#!/usr/bin/env python3
"""
Compare Ks distributions across different BLAST filtering thresholds.
Creates overlaid histograms and summary statistics.

Usage:
    python3 compare_ks_filtering.py input_dir output_prefix
    
Example:
    python3 analysis/Ks/compare_ks_filtering.py output/ks_results/ output/figures/ks_comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import re
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10

def extract_filter_params(filename):
    """Extract filtering parameters from filename."""
    # Pattern: ks_results_id50_qcov70_scov70.tsv
    match = re.search(r'id(\d+)_qcov(\d+)_scov(\d+)', filename)
    if match:
        return {
            'identity': int(match.group(1)),
            'qcov': int(match.group(2)),
            'scov': int(match.group(3)),
            'label': f"ID{match.group(1)}_Qcov{match.group(2)}_Scov{match.group(3)}"
        }
    elif 'filtered' in filename:
        return {
            'identity': 0,
            'qcov': 0,
            'scov': 0,
            'label': 'Baseline (filtered)'
        }
    else:
        return {
            'identity': 0,
            'qcov': 0,
            'scov': 0,
            'label': 'All pairs'
        }

def load_ks_files(input_dir):
    """Load all Ks results files and return as dict."""
    input_path = Path(input_dir)
    ks_files = list(input_path.glob("ks_results*.tsv"))
    
    if not ks_files:
        print(f"ERROR: No ks_results*.tsv files found in {input_dir}")
        sys.exit(1)
    
    data = {}
    for f in sorted(ks_files):
        params = extract_filter_params(f.name)
        df = pd.read_csv(f, sep='\t')
        # Filter to reasonable Ks range
        df_filtered = df[(df['ks'] > 0) & (df['ks'] <= 5)].copy()
        
        data[params['label']] = {
            'df': df_filtered,
            'params': params,
            'file': f.name
        }
        print(f"Loaded {f.name}: {len(df_filtered)} pairs (Ks 0-2)")
    
    return data

def fit_gmm_peaks(ks_values, n_components=2):
    """Fit Gaussian mixture model to detect peaks."""
    try:
        ks_array = ks_values.values.reshape(-1, 1)
        gmm = GaussianMixture(n_components=n_components, random_state=42, max_iter=1000)
        gmm.fit(ks_array)
        
        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        weights = gmm.weights_
        
        # Sort by mean
        sorted_idx = np.argsort(means)
        return {
            'means': means[sorted_idx],
            'stds': stds[sorted_idx],
            'weights': weights[sorted_idx]
        }
    except:
        return None

def plot_overlaid_distributions(data, output_prefix):
    """Create overlaid histogram plot."""
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Color palette
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(data)))
    
    # Plot each distribution
    for idx, (label, info) in enumerate(sorted(data.items(), 
                                                key=lambda x: x[1]['params']['identity'])):
        ks_vals = info['df']['ks']
        
        # Histogram with transparency
        ax.hist(ks_vals, bins=100, alpha=0.4, label=f"{label} (n={len(ks_vals)})",
                color=colors[idx], density=True, edgecolor='none')
        
        # Density line
        density = gaussian_kde(ks_vals)
        xs = np.linspace(0, 5, 300)
        ax.plot(xs, density(xs), color=colors[idx], linewidth=2, alpha=0.8)
    
    ax.set_xlabel('Ks (synonymous substitutions per site)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title('Ks Distribution Comparison Across Filtering Thresholds', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim(0, 5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_overlaid.png", dpi=300)
    plt.savefig(f"{output_prefix}_overlaid.pdf")
    print(f"✓ Saved overlaid plot: {output_prefix}_overlaid.png")
    plt.close()

def plot_subplots_grid(data, output_prefix):
    """Create subplot grid with individual distributions."""
    n_plots = len(data)
    n_cols = 2
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 5 * n_rows))
    axes = axes.flatten() if n_plots > 1 else [axes]
    
    for idx, (label, info) in enumerate(sorted(data.items(), 
                                                key=lambda x: x[1]['params']['identity'])):
        ax = axes[idx]
        ks_vals = info['df']['ks']
        
        # Histogram
        ax.hist(ks_vals, bins=100, alpha=0.7, color='steelblue', edgecolor='black', 
                linewidth=0.5, density=True)
        
        # Density line
        density = gaussian_kde(ks_vals)
        xs = np.linspace(0, 5, 300)
        ax.plot(xs, density(xs), 'r-', linewidth=2, label='Density')
        
        # Median line
        median_ks = ks_vals.median()
        ax.axvline(median_ks, color='darkgreen', linestyle='--', linewidth=2, 
                   label=f'Median: {median_ks:.3f}')
        
        # GMM peaks
        gmm_result = fit_gmm_peaks(ks_vals, n_components=2)
        if gmm_result:
            for peak_idx, mean in enumerate(gmm_result['means']):
                ax.axvline(mean, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
                ax.text(mean, ax.get_ylim()[1] * 0.9, f'{mean:.2f}', 
                       ha='center', fontsize=8, color='orange')
        
        ax.set_title(f"{label}\n(n={len(ks_vals)} pairs)", fontweight='bold')
        ax.set_xlabel('Ks')
        ax.set_ylabel('Density')
        ax.legend(fontsize=8)
        ax.set_xlim(0, 5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Hide unused subplots
    for idx in range(len(data), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_grid.png", dpi=300)
    plt.savefig(f"{output_prefix}_grid.pdf")
    print(f"✓ Saved grid plot: {output_prefix}_grid.png")
    plt.close()

def create_summary_table(data, output_prefix):
    """Create summary statistics table."""
    summary_data = []
    
    for label, info in sorted(data.items(), key=lambda x: x[1]['params']['identity']):
        ks_vals = info['df']['ks']
        
        # Fit GMM for peaks
        gmm_result = fit_gmm_peaks(ks_vals, n_components=2)
        
        row = {
            'Filter': label,
            'N_pairs': len(ks_vals),
            'Mean_Ks': ks_vals.mean(),
            'Median_Ks': ks_vals.median(),
            'SD_Ks': ks_vals.std(),
            'Min_Ks': ks_vals.min(),
            'Max_Ks': ks_vals.max(),
        }
        
        if gmm_result:
            row['Peak1_Ks'] = gmm_result['means'][0]
            row['Peak2_Ks'] = gmm_result['means'][1]
            row['Peak1_weight'] = gmm_result['weights'][0]
            row['Peak2_weight'] = gmm_result['weights'][1]
        else:
            row['Peak1_Ks'] = np.nan
            row['Peak2_Ks'] = np.nan
            row['Peak1_weight'] = np.nan
            row['Peak2_weight'] = np.nan
        
        summary_data.append(row)
    
    summary_df = pd.DataFrame(summary_data)
    
    # Save to TSV
    summary_df.to_csv(f"{output_prefix}_summary.tsv", sep='\t', index=False, float_format='%.4f')
    print(f"✓ Saved summary table: {output_prefix}_summary.tsv")
    
    # Create visual table
    fig, ax = plt.subplots(figsize=(16, len(summary_data) * 0.6 + 1))
    ax.axis('tight')
    ax.axis('off')
    
    # Format table
    table_data = summary_df.copy()
    for col in ['Mean_Ks', 'Median_Ks', 'SD_Ks', 'Peak1_Ks', 'Peak2_Ks', 
                'Peak1_weight', 'Peak2_weight']:
        if col in table_data.columns:
            table_data[col] = table_data[col].apply(lambda x: f'{x:.3f}' if pd.notna(x) else '-')
    table_data['N_pairs'] = table_data['N_pairs'].apply(lambda x: f'{x:,}')
    
    table = ax.table(cellText=table_data.values, colLabels=table_data.columns,
                    cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # Style header
    for i in range(len(table_data.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Alternate row colors
    for i in range(1, len(table_data) + 1):
        color = '#D9E1F2' if i % 2 == 0 else 'white'
        for j in range(len(table_data.columns)):
            table[(i, j)].set_facecolor(color)
    
    plt.title('Ks Distribution Summary Statistics by Filtering Threshold', 
              fontsize=14, fontweight='bold', pad=20)
    plt.savefig(f"{output_prefix}_summary_table.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved summary table image: {output_prefix}_summary_table.png")
    plt.close()
    
    return summary_df

def plot_combined_figure(data, output_prefix):
    """Create a comprehensive combined figure."""
    fig = plt.figure(figsize=(18, 10))
    
    # Main overlaid plot (top, large)
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2, rowspan=2)
    
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(data)))
    
    for idx, (label, info) in enumerate(sorted(data.items(), 
                                                key=lambda x: x[1]['params']['identity'])):
        ks_vals = info['df']['ks']
        ax1.hist(ks_vals, bins=100, alpha=0.35, label=f"{label} (n={len(ks_vals):,})",
                color=colors[idx], density=True, edgecolor='none')
        density = gaussian_kde(ks_vals)
        xs = np.linspace(0, 5, 300)
        ax1.plot(xs, density(xs), color=colors[idx], linewidth=2.5, alpha=0.9)
    
    ax1.set_xlabel('Ks (synonymous substitutions per site)', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.set_title('Ks Distribution Comparison Across Filtering Thresholds', 
                  fontsize=14, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax1.set_xlim(0, 5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.grid(alpha=0.3)
    
    # Statistics comparison plot (bottom left)
    ax2 = plt.subplot2grid((3, 2), (2, 0))
    
    labels_list = []
    means = []
    medians = []
    
    for label, info in sorted(data.items(), key=lambda x: x[1]['params']['identity']):
        labels_list.append(label.split('_')[0])  # Shortened label
        means.append(info['df']['ks'].mean())
        medians.append(info['df']['ks'].median())
    
    x = np.arange(len(labels_list))
    width = 0.35
    ax2.bar(x - width/2, means, width, label='Mean', alpha=0.8, color='steelblue')
    ax2.bar(x + width/2, medians, width, label='Median', alpha=0.8, color='coral')
    ax2.set_xlabel('Filter', fontsize=10)
    ax2.set_ylabel('Ks', fontsize=10)
    ax2.set_title('Mean vs Median Ks', fontsize=11, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels_list, rotation=45, ha='right', fontsize=8)
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3, axis='y')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # Pair count comparison (bottom right)
    ax3 = plt.subplot2grid((3, 2), (2, 1))
    
    counts = [len(info['df']) for label, info in sorted(data.items(), 
                                                         key=lambda x: x[1]['params']['identity'])]
    ax3.bar(labels_list, counts, alpha=0.8, color='darkgreen')
    ax3.set_xlabel('Filter', fontsize=10)
    ax3.set_ylabel('Number of pairs', fontsize=10)
    ax3.set_title('Gene Pair Count by Filter', fontsize=11, fontweight='bold')
    ax3.set_xticklabels(labels_list, rotation=45, ha='right', fontsize=8)
    ax3.grid(alpha=0.3, axis='y')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    # Format y-axis with comma separator
    ax3.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_combined.png", dpi=300)
    plt.savefig(f"{output_prefix}_combined.pdf")
    print(f"✓ Saved combined figure: {output_prefix}_combined.png")
    plt.close()

def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_prefix = sys.argv[2]
    
    print("=== Ks Distribution Comparison ===\n")
    
    # Load all Ks files
    data = load_ks_files(input_dir)
    print(f"\nLoaded {len(data)} datasets")
    
    # Create output directory
    output_path = Path(output_prefix).parent
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    print("\nGenerating visualizations...")
    plot_overlaid_distributions(data, output_prefix)
    plot_subplots_grid(data, output_prefix)
    plot_combined_figure(data, output_prefix)
    
    # Generate summary table
    print("\nGenerating summary statistics...")
    summary_df = create_summary_table(data, output_prefix)
    
    print("\n=== Summary ===")
    print(summary_df[['Filter', 'N_pairs', 'Mean_Ks', 'Median_Ks', 'Peak1_Ks', 'Peak2_Ks']].to_string(index=False))
    
    print("\n✓ All visualizations complete!")
    print(f"  Output files: {output_prefix}_*")

if __name__ == "__main__":
    main()
