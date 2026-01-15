#!/usr/bin/env python3
"""
Plot duplication age and Ks distribution for our generated anchors.
Similar to MCScanX anchor plot but using our collinearity-detected anchors.

Input: Anchors TSV file (gene1, gene2, ks columns)
Output: Age/Ks distribution plots
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from sklearn.mixture import GaussianMixture
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Anchors TSV file (gene1, gene2, ks)')
parser.add_argument('-o', '--output', required=True, help='Output prefix for plots')
parser.add_argument('--stats-dir', default=None, help='Directory for statistics output')
args = parser.parse_args()

output_prefix = args.output
stats_dir = Path(args.stats_dir) if args.stats_dir else Path(output_prefix).parent / 'statistics'
stats_dir.mkdir(parents=True, exist_ok=True)

# Load anchors
print(f"\nLoading anchors from: {args.input}")
df = pd.read_csv(args.input, sep='\t')
print(f"  Total anchor pairs: {len(df)}")

# Filter Ks (0 < Ks <= 2)
df = df[df['ks'].between(0, 2, inclusive='both') & df['ks'].notna()].copy()
print(f"  After Ks filter (0 < Ks ≤ 2): {len(df)}")

if len(df) == 0:
    print("No pairs remaining after filtering. Exiting.")
    exit(1)

# Calculate age using lambda
lambda_rate = 6.1e-9  # substitutions per site per year (legumes)
df['Age_MYA'] = df['ks'] / (2 * lambda_rate) / 1e6

# === Statistics ===
print("\n=== Ks Statistics ===")
print(f"Mean Ks: {df['ks'].mean():.4f}")
print(f"Median Ks: {df['ks'].median():.4f}")
print(f"SD Ks: {df['ks'].std():.4f}")
print(f"Min Ks: {df['ks'].min():.4f}")
print(f"Max Ks: {df['ks'].max():.4f}")

print("\n=== Age Statistics (MYA) ===")
print(f"Mean Age: {df['Age_MYA'].mean():.1f}")
print(f"Median Age: {df['Age_MYA'].median():.1f}")
print(f"SD Age: {df['Age_MYA'].std():.1f}")

# Save summary
stats_summary = pd.DataFrame({
    'Metric': ['Total_Pairs', 'Mean_Ks', 'Median_Ks', 'Std_Ks', 'Min_Ks', 'Max_Ks',
               'Mean_Age_MYA', 'Median_Age_MYA', 'Std_Age_MYA', 'Min_Age_MYA', 'Max_Age_MYA'],
    'Value': [len(df), df['ks'].mean(), df['ks'].median(), df['ks'].std(), df['ks'].min(), df['ks'].max(),
              df['Age_MYA'].mean(), df['Age_MYA'].median(), df['Age_MYA'].std(), df['Age_MYA'].min(), df['Age_MYA'].max()]
})
stats_summary.to_csv(stats_dir / 'anchors_summary_stats.tsv', sep='\t', index=False)
print(f"\nStatistics saved to: {stats_dir / 'anchors_summary_stats.tsv'}")

# === WGD Peak Detection ===
print("\n=== WGD Peak Detection ===")
try:
    ks_wgd = df[df['ks'] < 2.5]['ks'].values.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=42, max_iter=1000)
    gmm.fit(ks_wgd)
    
    means = gmm.means_.flatten()
    variances = gmm.covariances_.flatten()
    weights = gmm.weights_
    
    sorted_idx = np.argsort(means)
    means = means[sorted_idx]
    variances = variances[sorted_idx]
    stds = np.sqrt(variances)
    weights = weights[sorted_idx]
    
    print("\nMixture model results:")
    print(f"Component 1: μ = {means[0]:.4f}, σ = {stds[0]:.4f}, weight = {weights[0]:.3f}")
    print(f"Component 2: μ = {means[1]:.4f}, σ = {stds[1]:.4f}, weight = {weights[1]:.3f}")
    
    age1 = means[0] / (2 * lambda_rate) / 1e6
    age2 = means[1] / (2 * lambda_rate) / 1e6
    
    print(f"\nEstimated WGD ages:")
    print(f"  Peak 1: {age1:.1f} MYA (Ks = {means[0]:.4f})")
    print(f"  Peak 2: {age2:.1f} MYA (Ks = {means[1]:.4f})")
    
    peak_info = pd.DataFrame({
        'Peak': [1, 2],
        'Ks_mean': means,
        'Ks_sd': stds,
        'Proportion': weights,
        'Age_MYA': [age1, age2]
    })
    peak_info.to_csv(stats_dir / 'anchors_ks_wgd_peaks.tsv', sep='\t', index=False)
    print(f"\nPeaks saved to: {stats_dir / 'anchors_ks_wgd_peaks.tsv'}")
    
except Exception as e:
    print(f"Warning: Peak detection failed: {str(e)}")

# === Plot age distribution ===
print("\nGenerating age distribution plot...")

fig, ax = plt.subplots(figsize=(10, 6))

# Histogram
ax.hist(df['Age_MYA'], bins=100, alpha=0.7, color='#084594', 
        label=f'Our Anchors (N={len(df):,})', density=True)

# Density curve
density = gaussian_kde(df['Age_MYA'])
xs = np.linspace(df['Age_MYA'].min(), df['Age_MYA'].max(), 200)
ax.plot(xs, density(xs), color='#7A0177', linewidth=2, label='Density')

# Known WGD events
ax.axvline(13, color='#7A0177', linestyle='--', linewidth=2)
ax.axvline(59, color='#7A0177', linestyle='--', linewidth=2)
ax.text(13, ax.get_ylim()[1] * 0.95, '~13 MYA\n(known WGD)', 
        ha='left', va='top', color='#7A0177', fontsize=9)
ax.text(59, ax.get_ylim()[1] * 0.95, '~59 MYA\n(known WGD)', 
        ha='left', va='top', color='#7A0177', fontsize=9)

ax.set_xlabel('Age (Million Years Ago)')
ax.set_ylabel('Density')
ax.set_title('Duplication Age Distribution (Our Detected Anchors)', fontweight='bold', pad=20)
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Secondary x-axis for Ks
def age_to_ks(age):
    return age * 2 * lambda_rate * 1e6
def ks_to_age(ks):
    return ks / (2 * lambda_rate) / 1e6

secax = ax.secondary_xaxis('top', functions=(age_to_ks, ks_to_age))
secax.set_xlabel('Ks (synonymous substitutions per site)')
ks_max = df['ks'].max()
secax.set_xticks(np.arange(0, ks_max + 0.1, 0.1))

plt.tight_layout()
plt.savefig(f'{output_prefix}_age_ks_plot.png', dpi=300)
plt.savefig(f'{output_prefix}_age_ks_plot.pdf')
plt.close()

print(f"\nPlots saved:")
print(f"  {output_prefix}_age_ks_plot.png")
print(f"  {output_prefix}_age_ks_plot.pdf")
print("\n=== Complete ===")
