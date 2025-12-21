#!/usr/bin/env python3
"""
Plot duplication age and Ks distribution for MCScanX-anchored gene pairs.
Inputs:
  --in  Path to anchors_with_ks.tsv (from MCScanX annotation)
  --out Output prefix for plots (PNG/PDF)
  --stats-dir Directory for statistics output (TSV files)
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy import stats
from sklearn.mixture import GaussianMixture
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, help='Anchors with Ks TSV file')
parser.add_argument('-o', required=True, help='Output prefix for plots')
parser.add_argument('--stats-dir', default=None, help='Directory for statistics output (default: same dir as plots)')
args = parser.parse_args()

input_file = args.i
output_prefix = args.o
stats_dir = Path(args.stats_dir) if args.stats_dir else Path(output_prefix).parent / 'statistics'
stats_dir.mkdir(parents=True, exist_ok=True)

# Load anchored Ks data
df = pd.read_csv(input_file, sep='\t')
df = df[df['ks'].between(0, 2, inclusive='both') & df['ks'].notna()]

# Calculate age
lambda_rate = 6.1e-9
if 'Age_MYA' not in df.columns:
    df['Age_MYA'] = df['ks'] / (2 * lambda_rate) / 1e6

# === Basic Statistics ===
print("\n=== Anchored Pairs - Ks Statistics ===")
print(df['ks'].describe())
print(f"Mean Ks: {df['ks'].mean():.3f}")
print(f"Median Ks: {df['ks'].median():.3f}")
print(f"SD Ks: {df['ks'].std():.3f}")

print("\n=== Anchored Pairs - Age Statistics ===")
print(df['Age_MYA'].describe())
print(f"Mean Age: {df['Age_MYA'].mean():.1f} MYA")
print(f"Median Age: {df['Age_MYA'].median():.1f} MYA")
print(f"SD Age: {df['Age_MYA'].std():.1f} MYA")

# Save summary statistics
stats_summary = pd.DataFrame({
    'Metric': ['Total_Pairs', 'Mean_Ks', 'Median_Ks', 'Std_Ks', 'Min_Ks', 'Max_Ks',
               'Mean_Age_MYA', 'Median_Age_MYA', 'Std_Age_MYA', 'Min_Age_MYA', 'Max_Age_MYA'],
    'Value': [len(df), df['ks'].mean(), df['ks'].median(), df['ks'].std(), df['ks'].min(), df['ks'].max(),
              df['Age_MYA'].mean(), df['Age_MYA'].median(), df['Age_MYA'].std(), df['Age_MYA'].min(), df['Age_MYA'].max()]
})
stats_summary.to_csv(stats_dir / 'anchored_pairs_summary_stats.tsv', sep='\t', index=False)

# Save filtered data
df.to_csv(stats_dir / 'anchored_pairs_filtered.tsv', sep='\t', index=False)

print(f"\nStatistics saved to: {stats_dir}")

# === WGD Peak Detection ===
print("\n=== WGD Peak Detection ===")

try:
    # Focus on Ks < 1.4 for WGD detection
    ks_wgd = df[df['ks'] < 2.5]['ks'].values.reshape(-1, 1)
    
    # Fit 2-component mixture
    gmm = GaussianMixture(n_components=2, random_state=42, max_iter=1000)
    gmm.fit(ks_wgd)
    
    # Get parameters
    means = gmm.means_.flatten()
    variances = gmm.covariances_.flatten()
    weights = gmm.weights_
    
    # Sort by mean
    sorted_idx = np.argsort(means)
    means = means[sorted_idx]
    variances = variances[sorted_idx]
    stds = np.sqrt(variances)
    weights = weights[sorted_idx]
    
    print("\nMixture model results:")
    print(f"Component 1: μ = {means[0]:.3f}, σ = {stds[0]:.3f}, λ = {weights[0]:.3f}")
    print(f"Component 2: μ = {means[1]:.3f}, σ = {stds[1]:.3f}, λ = {weights[1]:.3f}")
    
    # Convert Ks to age (MYA)
    lambda_rate_legume = 6.1e-9  # substitutions per site per year (legumes)
    age1 = means[0] / (2 * lambda_rate_legume) / 1e6
    age2 = means[1] / (2 * lambda_rate_legume) / 1e6
    
    print(f"\nEstimated ages (λ = {lambda_rate_legume:.1e}):")
    print(f"  Peak 1: {age1:.1f} MYA")
    print(f"  Peak 2: {age2:.1f} MYA")
    
    # Save peak information
    peak_info = pd.DataFrame({
        'Peak': [1, 2],
        'Ks_mean': means,
        'Ks_sd': stds,
        'Proportion': weights,
        'Age_MYA': [age1, age2]
    })
    
    peak_info.to_csv(stats_dir / 'ks_wgd_peaks.tsv', sep='\t', index=False)
    print(f"WGD peaks saved to: {stats_dir / 'ks_wgd_peaks.tsv'}")
    
except Exception as e:
    print("Warning: Mixture model fitting failed. Skipping peak detection.")
    print(f"Error: {str(e)}")

# Plot age histogram with secondary Ks axis
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(df['Age_MYA'], bins=100, alpha=0.7, color='#084594', label='Anchored Data', density=True)
# Add density curve (KDE)
density = gaussian_kde(df['Age_MYA'])
xs = np.linspace(df['Age_MYA'].min(), df['Age_MYA'].max(), 200)
ax.plot(xs, density(xs), color='#7A0177', linewidth=2, label='Density')

# Highlight known WGD events
ax.axvline(13, color='#7A0177', linestyle='--', linewidth=2)
ax.axvline(59, color='#7A0177', linestyle='--', linewidth=2)
ax.text(13, ax.get_ylim()[1] * 0.95, '~13 MYA\n(known WGD)', 
    ha='left', va='top', color='#7A0177', fontsize=9)
ax.text(59, ax.get_ylim()[1] * 0.95, '~59 MYA\n(known WGD)', 
    ha='left', va='top', color='#7A0177', fontsize=9)

ax.set_xlabel('Age (Million Years Ago)')
ax.set_ylabel('Count')
ax.set_title('Duplication Age Distribution (Anchored Pairs)', fontweight='bold', pad=20)
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

def age_to_ks(age):
    return age * 2 * lambda_rate * 1e6
def ks_to_age(ks):
    return ks / (2 * lambda_rate) / 1e6

secax = ax.secondary_xaxis('top', functions=(age_to_ks, ks_to_age))
secax.set_xlabel('Ks (synonymous substitutions per site)')
ks_max = df['ks'].max()
secax.set_ticks(np.arange(0, ks_max+0.1, 0.1))

plt.tight_layout()
plt.savefig(f'{output_prefix}_age_ks_plot.png', dpi=300)
plt.savefig(f'{output_prefix}_age_ks_plot.pdf')
plt.close()
