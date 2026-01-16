import os
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind, chi2_contingency
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def cohens_d(x, y):
    """Calculate Cohen's d effect size"""
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx-1)*np.var(x, ddof=1) + (ny-1)*np.var(y, ddof=1)) / dof)
    return (np.mean(x) - np.mean(y)) / pooled_std

def main():
    parser = argparse.ArgumentParser(description='Analyze TAG age distribution using Ks values')
    parser.add_argument('--tag-pairs', default='analysis/duplication_types/TAGs/TAGs_distance_100000bp.tsv',
                       help='TAG pairs file from identify_TAGs.py')
    parser.add_argument('--ks-file', default='output/pipeline2/ks_results/ks_results.tsv',
                       help='Ks results file')
    parser.add_argument('--output-dir', default='analysis/duplication_types/TAGs/plots',
                       help='Output directory for plots')
    parser.add_argument('--stats-dir', default='analysis/duplication_types/TAGs/stats',
                       help='Output directory for statistics')
    
    args = parser.parse_args()
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.stats_dir, exist_ok=True)
    
    # === Data Loading ===
    print(f"Loading TAG pairs from {args.tag_pairs}...")
    tag_pairs = pd.read_csv(args.tag_pairs, sep='\t')
    print(f"  Loaded {len(tag_pairs)} TAG pairs")
    
    print(f"Loading Ks values from {args.ks_file}...")
    ks = pd.read_csv(args.ks_file, sep='\t')
    print(f"  Loaded {len(ks)} Ks values")
    
    # Create a set of TAG gene pairs for faster lookup
    tag_gene_pairs = set()
    for _, row in tag_pairs.iterrows():
        # Add both orderings of the pair
        tag_gene_pairs.add((row['Gene1'], row['Gene2']))
        tag_gene_pairs.add((row['Gene2'], row['Gene1']))
    
    # Mark Ks values as TAG or Non-TAG
    # Assuming ks file has columns like 'gene1', 'gene2', 'dS' or similar
    # Adjust column names as needed
    ks_cols = ks.columns.tolist()
    print(f"Ks columns: {ks_cols}")
    gene1_col = 'gene1' if 'gene1' in ks_cols else 'Gene1'
    gene2_col = 'gene2' if 'gene2' in ks_cols else 'Gene2'
    
    ks_col = 'ks'
    
    ks['Is_TAG'] = ks.apply(lambda row: (row[gene1_col], row[gene2_col]) in tag_gene_pairs, axis=1)
    ks['Duplication_Type'] = np.where(ks['Is_TAG'], 'TAG', 'Non-TAG')
    
    print(f"  TAG pairs in Ks data: {ks['Is_TAG'].sum()}")
    print(f"  Non-TAG pairs: {(~ks['Is_TAG']).sum()}")
    
    # Filter valid Ks values
    ks_valid = ks[ks[ks_col].notna() & (ks[ks_col] > 0) & (ks[ks_col] < 2)]
    
    # === Statistical Tests ===
    ks_tag = ks_valid[ks_valid['Duplication_Type'] == 'TAG'][ks_col].dropna()
    ks_nontag = ks_valid[ks_valid['Duplication_Type'] == 'Non-TAG'][ks_col].dropna()
    
    print(f"\n=== Statistical Analysis ===")
    print(f"TAG Ks values: {len(ks_tag)} (mean: {ks_tag.mean():.4f}, median: {ks_tag.median():.4f})")
    print(f"Non-TAG Ks values: {len(ks_nontag)} (mean: {ks_nontag.mean():.4f}, median: {ks_nontag.median():.4f})")
    
    # Mann-Whitney U (Wilcoxon rank-sum)
    wilcox_stat, wilcox_p = mannwhitneyu(ks_tag, ks_nontag, alternative='less')
    # T-test
    t_stat, t_p = ttest_ind(ks_tag, ks_nontag, equal_var=False)
    # Cohen's d
    cohen_d = cohens_d(ks_tag, ks_nontag)
    
    print(f"\nMann-Whitney U test: statistic={wilcox_stat:.2f}, p-value={wilcox_p:.2e}")
    print(f"T-test: statistic={t_stat:.4f}, p-value={t_p:.2e}")
    print(f"Cohen's d: {cohen_d:.4f}")
    
    # === Density Plot ===
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.kdeplot(ks_tag, label='TAG', fill=True, color="#E74C3C", alpha=0.6)
    sns.kdeplot(ks_nontag, label='Non-TAG', fill=True, color="#2980B9", alpha=0.6)
    plt.axvline(np.median(ks_tag), color="#E74C3C", linestyle='--', label='TAG Median')
    plt.axvline(np.median(ks_nontag), color="#2980B9", linestyle='--', label='Non-TAG Median')
    plt.title("Ks Distribution: TAG vs Non-TAG", fontsize=16, fontweight='bold')
    plt.xlabel("Ks (synonymous substitutions per site)")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "ks_TAG_vs_nonTAG_density.pdf"))
    plt.savefig(os.path.join(args.output_dir, "ks_TAG_vs_nonTAG_density.png"), dpi=300)
    plt.close()
    
    # === Orientation Analysis ===
    if 'Same_Orientation' in tag_pairs.columns:
        print(f"\n=== Orientation Analysis ===")
        orientation_counts = tag_pairs['Same_Orientation'].value_counts()
        same = orientation_counts.get(True, 0)
        total = orientation_counts.sum()
        print(f"Same orientation: {same} / {total} ({same/total*100:.1f}%)")
        
        # Chi-square test
        obs = [orientation_counts.get(True, 0), orientation_counts.get(False, 0)]
        chi2, p_chi2, _, _ = chi2_contingency([obs, [total/2, total/2]])
        print(f"Chi-square test: χ²={chi2:.2f}, p-value={p_chi2:.2e}")
        
        # Bar plot
        plt.figure(figsize=(8, 6))
        tag_pairs['Same_Orientation_str'] = tag_pairs['Same_Orientation'].astype(str)
        sns.countplot(x='Same_Orientation_str', data=tag_pairs, hue='Same_Orientation_str', legend=False,
                      palette={'False': "#95A5A6", 'True': "#27AE60"})
        plt.title(f"TAG Gene Orientation\nχ² test: p = {p_chi2:.2e}", fontsize=14, fontweight='bold')
        plt.xlabel("Orientation")
        plt.ylabel("Number of TAG pairs")
        plt.xticks([0, 1], ['Opposite', 'Same'])
        plt.tight_layout()
        plt.savefig(os.path.join(args.output_dir, "TAG_orientation.pdf"))
        plt.savefig(os.path.join(args.output_dir, "TAG_orientation.png"), dpi=300)
        plt.close()
    
    # === Save Results ===
    print(f"\n=== Saving Results ===")
    
    # Statistical test results
    test_results = pd.DataFrame({
        'Test': ["Mann-Whitney U", "T-test", "Cohen's d"],
        'Statistic': [wilcox_stat, t_stat, cohen_d],
        'P_value': [wilcox_p, t_p, np.nan],
        'Conclusion': [
            "TAGs younger" if wilcox_p < 0.05 else "No difference",
            "TAGs younger" if t_p < 0.05 else "No difference",
            "Medium/Large effect" if abs(cohen_d) > 0.5 else "Small effect"
        ]
    })
    test_results.to_csv(os.path.join(args.stats_dir, "TAG_age_tests.tsv"), sep='\t', index=False)
    
    # Summary statistics
    summary_stats = ks_valid.groupby('Duplication_Type')[ks_col].agg([
        ('N', 'count'),
        ('Mean_Ks', 'mean'),
        ('Median_Ks', 'median'),
        ('SD_Ks', 'std'),
        ('Q25', lambda x: x.quantile(0.25)),
        ('Q75', lambda x: x.quantile(0.75))
    ]).reset_index()
    summary_stats.to_csv(os.path.join(args.stats_dir, "TAG_ks_summary.tsv"), sep='\t', index=False)
    
    print(f"Plots saved to: {args.output_dir}")
    print(f"Statistics saved to: {args.stats_dir}")
    print("\n=== Analysis Complete ===")

if __name__ == '__main__':
    main()
