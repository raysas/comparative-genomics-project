#!/usr/bin/env python3
"""
Ks distribution analysis and WGD detection
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.mixture import GaussianMixture
import os
import sys

# Configuration
KS_FILE = sys.argv[1] if len(sys.argv) > 1 else "../../output/ks_results/ks_results_filtered.tsv"
OUTPUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "figures/"
STATS_DIR = sys.argv[3] if len(sys.argv) > 3 else "statistics/"

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(STATS_DIR, exist_ok=True)

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 10

# Load Ks data
print("Loading Ks data...")
ks_data = pd.read_csv(KS_FILE, sep='\t')
print(f"  Loaded {len(ks_data)} gene pairs")

# Quality filtering
print("\nApplying quality filters...")
print(f"  Original pairs: {len(ks_data)}")

# Filter saturated pairs (Ks > 5)
ks_filtered = ks_data[(ks_data['ks'] <=1.4) & (ks_data['ks'] > 0)].copy()
print(f"  After Ks filter (0 < Ks <= 5): {len(ks_filtered)}")

# Additional QC filters (adjust column names as needed)
# ks_filtered = ks_filtered[(ks_filtered['Ka_Ks'] < 5) & (ks_filtered['Alignment_Length'] > 100)]

# === Basic Statistics ===
print("\n=== Ks Statistics ===")
print(ks_filtered['ks'].describe())
print(f"Mean Ks: {ks_filtered['ks'].mean():.3f}")
print(f"Median Ks: {ks_filtered['ks'].median():.3f}")
print(f"SD Ks: {ks_filtered['ks'].std():.3f}")

# === Ks Distribution Plot ===
print("\nGenerating Ks distribution plot...")

fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(ks_filtered['ks'], bins=100, density=True, alpha=0.7, color='steelblue', label='Data')

# Add density line
from scipy.stats import gaussian_kde
density = gaussian_kde(ks_filtered['ks'])
xs = np.linspace(ks_filtered['ks'].min(), ks_filtered['ks'].max(), 200)
ax.plot(xs, density(xs), 'r-', linewidth=2, label='Density')

# Add median line
median_ks = ks_filtered['ks'].median()
ax.axvline(median_ks, color='darkgreen', linestyle='--', linewidth=2, label=f'Median: {median_ks:.3f}')

ax.set_xlabel('Ks (synonymous substitutions per site)')
ax.set_ylabel('Density')
ax.set_title(f'Ks Distribution - All Duplicated Genes\nN = {len(ks_filtered)} gene pairs', fontweight='bold')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'ks_distribution_all.pdf'))
plt.savefig(os.path.join(OUTPUT_DIR, 'ks_distribution_all.png'), dpi=300)
plt.close()

# === WGD Peak Detection ===
print("\n=== WGD Peak Detection ===")

# Expected WGD events in Glycine max:
# - ~13 MYA: Ks ~ 0.1-0.2 (λ = 6.5e-9, Age = Ks / (2*λ) )
# - ~59 MYA: Ks ~ 0.5-0.6

# Fit mixture of normal distributions
print("Fitting Gaussian mixture model...")

try:
    # Focus on Ks < 2 for WGD detection
    ks_wgd = ks_filtered[ks_filtered['ks'] < 1.4]['ks'].values.reshape(-1, 1)
    
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
    
    # Plot with mixture components
    x_vals = np.linspace(0, 2, 200)
    
    # Calculate mixture density
    comp1 = weights[0] * stats.norm.pdf(x_vals, means[0], stds[0])
    comp2 = weights[1] * stats.norm.pdf(x_vals, means[1], stds[1])
    mixture_density = comp1 + comp2
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(ks_wgd, bins=50, density=True, alpha=0.5, color='steelblue', label='Data')
    ax.plot(x_vals, mixture_density, 'r-', linewidth=2, label='Mixture model')
    ax.plot(x_vals, comp1, 'orange', linestyle='--', linewidth=1.5, label=f'Component 1 (μ={means[0]:.3f})')
    ax.plot(x_vals, comp2, 'purple', linestyle='--', linewidth=1.5, label=f'Component 2 (μ={means[1]:.3f})')
    ax.axvline(means[0], color='orange', linestyle=':', linewidth=1)
    ax.axvline(means[1], color='purple', linestyle=':', linewidth=1)
    
    ax.set_xlabel('Ks (synonymous substitutions per site)')
    ax.set_ylabel('Density')
    ax.set_title('Ks Distribution with WGD Peak Detection\nGaussian Mixture Model (2 components)', fontweight='bold')
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'ks_wgd_peaks.pdf'))
    plt.savefig(os.path.join(OUTPUT_DIR, 'ks_wgd_peaks.png'), dpi=300)
    plt.close()
    
    # Convert Ks to age (MYA)
    lambda_rate = 6.5e-9  # substitutions per site per year (legumes)
    age1 = means[0] / (2 * lambda_rate) / 1e6  # Convert to millions of years
    age2 = means[1] / (2 * lambda_rate) / 1e6
    
    print(f"\nEstimated ages (λ = {lambda_rate:.1e}):")
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
    
    peak_info.to_csv(os.path.join(STATS_DIR, 'ks_wgd_peaks.tsv'), sep='\t', index=False)
    
except Exception as e:
    print("Warning: Mixture model fitting failed. Skipping peak detection.")
    print(f"Error: {str(e)}")

# === Age Conversion ===
print("\nConverting Ks to divergence time...")

lambda_rate = 4.1e-9  # substitutions per site per year
ks_filtered['Age_MYA'] = ks_filtered['ks'] / (2 * lambda_rate) / 1e6

# Age distribution plot
age_data = ks_filtered[ks_filtered['Age_MYA'] < 200]

fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(age_data['Age_MYA'], bins=100, density=True, alpha=0.7, color='darkgreen', label='Data')

# Add density line
density = gaussian_kde(age_data['Age_MYA'])
xs = np.linspace(age_data['Age_MYA'].min(), age_data['Age_MYA'].max(), 200)
ax.plot(xs, density(xs), 'r-', linewidth=2, label='Density')

# Add known WGD events
ax.axvline(13, color='blue', linestyle='--', linewidth=2)
ax.axvline(59, color='blue', linestyle='--', linewidth=2)
ax.text(13, ax.get_ylim()[1] * 0.95, '~13 MYA\n(known WGD)', 
        ha='left', va='top', color='blue', fontsize=9)
ax.text(59, ax.get_ylim()[1] * 0.95, '~59 MYA\n(known WGD)', 
        ha='left', va='top', color='blue', fontsize=9)

ax.set_xlabel('Age (Million Years Ago)')
ax.set_ylabel('Density')
ax.set_title(f'Duplication Age Distribution\nConversion: Age (MYA) = Ks / (2 × λ), λ = {lambda_rate:.1e}', 
             fontweight='bold')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'age_distribution.pdf'))
plt.savefig(os.path.join(OUTPUT_DIR, 'age_distribution.png'), dpi=300)
plt.close()

# === Save filtered Ks data ===
print("\nSaving filtered Ks data...")
ks_filtered.to_csv(os.path.join(STATS_DIR, 'ks_filtered.tsv'), sep='\t', index=False)

print("\n=== Analysis Complete ===")
print(f"Plots saved to: {OUTPUT_DIR}")
print(f"Statistics saved to: {STATS_DIR}")
