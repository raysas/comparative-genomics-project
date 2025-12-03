#!/usr/bin/env python3
"""
Plot Ks distribution to detect whole genome duplication (WGD) events
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_ks_wgd(tsv_file, output_dir, max_ks=2.0, sample_size=100000):
    """Create Ks distribution plot optimized for WGD detection"""
    
    print(f"Loading {tsv_file}...")
    df = pd.read_csv(tsv_file, sep='\t')
    
    # Filter
    df = df[(df['ks'] >= 0) & (df['ks'] <= max_ks)].copy()
    
    # Sample for speed if needed
    # if len(df) > sample_size:
    #     df = df.sample(n=sample_size, random_state=42)
    #     print(f"Sampled {sample_size} records (total: {len(df)})")
    # else:
    #     print(f"Using all {len(df)} records")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # High-resolution histogram for WGD detection
    counts, bins, patches = ax.hist(df['ks'], bins=200, edgecolor='black', 
                                     alpha=0.7, color='steelblue', linewidth=0.5)
    
    # Add density overlay
    from scipy import stats
    density = stats.gaussian_kde(df['ks'])
    xs = np.linspace(df['ks'].min(), df['ks'].max(), 200)
    ax2 = ax.twinx()
    ax2.plot(xs, density(xs), 'r-', linewidth=2.5, label='Density estimate')
    ax2.set_ylabel('Density', fontsize=12, fontweight='bold', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    
    ax.set_xlabel('Ks (Synonymous Substitution Rate)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=13, fontweight='bold')
    ax.set_title('Ks Distribution - WGD Detection\n(Multiple peaks indicate duplication events)', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add statistics box
    mean_ks = df['ks'].mean()
    median_ks = df['ks'].median()
    textstr = f'Mean Ks: {mean_ks:.3f}\nMedian Ks: {median_ks:.3f}\nN: {len(df):,}'
    ax.text(0.98, 0.97, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    output_file = Path(output_dir) / 'ks_distribution_wgd.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()
    
    # Print statistics for WGD analysis
    print("\n" + "="*60)
    print("Ks STATISTICS FOR WGD DETECTION")
    print("="*60)
    print(f"Mean Ks:      {mean_ks:.4f}")
    print(f"Median Ks:    {median_ks:.4f}")
    print(f"Std Dev:      {df['ks'].std():.4f}")
    print(f"Min Ks:       {df['ks'].min():.4f}")
    print(f"Max Ks:       {df['ks'].max():.4f}")
    print(f"\nPercentiles:")
    for q in [0.05, 0.25, 0.5, 0.75, 0.95]:
        print(f"  {q*100:>2.0f}th: {df['ks'].quantile(q):.4f}")
    
    # Find peaks using simple method
    hist, bin_edges = np.histogram(df['ks'], bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Find local maxima
    peaks = []
    for i in range(1, len(hist)-1):
        if hist[i] > hist[i-1] and hist[i] > hist[i+1]:
            peaks.append((bin_centers[i], hist[i]))
    
    if peaks:
        peaks.sort(key=lambda x: x[1], reverse=True)
        print(f"\nDetected {len(peaks)} potential WGD peaks:")
        for i, (ks_val, count) in enumerate(peaks[:5], 1):
            print(f"  Peak {i}: Ks ≈ {ks_val:.3f} (height: {count:.0f})")
    
    print("="*60)
    print("\nINTERPRETATION:")
    print("- Single sharp peak: Recent single duplication event")
    print("- Multiple peaks: Multiple duplication events at different times")
    print("- Broad distribution: Continuous duplication or ancient event")
    print("="*60)

if __name__ == '__main__':
    # take input file from command line argument
    import sys
    tsv_file = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(exist_ok=True)
    
    plot_ks_wgd(str(tsv_file), output_dir, max_ks=1.5, sample_size=100000)
