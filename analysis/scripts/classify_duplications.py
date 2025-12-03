#!/usr/bin/env python3
"""
Classify gene duplications into different types:
- TAGs (Tandem Arrayed Genes): Same chromosome, ≤10 genes or ≤100kb apart
- Proximal: Same chromosome, >100kb but <500kb apart  
- Dispersed: Different chromosomes OR >500kb apart
- WGD: Based on Ks peak ranges (~0.1-0.2 for 13MYA, ~0.8-1.0 for 59MYA)

Usage:
    python3 classify_duplications.py
    
Outputs:
    - gene_pairs_classified.tsv: All gene pairs with duplication type annotation
    - duplication_type_summary.tsv: Summary statistics by type
    - duplication_type_plots.png: Visualization of classification
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
KS_FILE = "output/ks_results/ks_results_filtered.tsv"
if sys.argv and len(sys.argv) > 1:
    KS_FILE = sys.argv[1]
OUTPUT_DIR = "output/duplication_classification/"

# Classification parameters
TAG_GENE_DISTANCE = 10  # Maximum genes apart for TAGs
TAG_BP_DISTANCE = 100000  # 100kb
PROXIMAL_BP_MIN = 100000  # 100kb
PROXIMAL_BP_MAX = 500000  # 500kb
WGD_PEAK1_MIN = 0.1  # 13 MYA event
WGD_PEAK1_MAX = 0.2
WGD_PEAK2_MIN = 0.45  # 59 MYA event
WGD_PEAK2_MAX = 0.55

# Create output directory
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

def load_data():
    """Load protein info and Ks results."""
    print("Loading data...")
    
    # Load protein/gene information
    protein_df = pd.read_csv(PROTEIN_INFO)
    print(f"  ✓ Loaded {len(protein_df)} genes")
    
    # Create gene position mapping
    gene_info = protein_df.set_index('peptide_id')[['chromosome', 'start_pos', 'end_pos']].to_dict('index')
    
    # Load Ks results
    ks_df = pd.read_csv(KS_FILE, sep='\t')
    print(f"  ✓ Loaded {len(ks_df)} gene pairs with Ks")
    
    return gene_info, ks_df

def calculate_gene_order(protein_df):
    """Calculate gene order on each chromosome."""
    print("Calculating gene order on chromosomes...")
    
    gene_order = {}
    for chrom in protein_df['chromosome'].unique():
        chr_genes = protein_df[protein_df['chromosome'] == chrom].sort_values('start_pos')
        for idx, gene_id in enumerate(chr_genes['peptide_id']):
            gene_order[gene_id] = {
                'chromosome': chrom,
                'order': idx,
                'position': chr_genes[chr_genes['peptide_id'] == gene_id]['start_pos'].values[0]
            }
    
    print(f"  ✓ Calculated order for {len(gene_order)} genes")
    return gene_order

def classify_duplication_type(row, gene_info, gene_order):
    """Classify a gene pair into duplication type."""
    gene1 = row['gene1']
    gene2 = row['gene2']
    ks = row['ks']
    
    # Get gene info
    if gene1 not in gene_info or gene2 not in gene_info:
        return 'Unknown', np.nan, np.nan
    
    chr1 = gene_info[gene1]['chromosome']
    chr2 = gene_info[gene2]['chromosome']
    pos1 = gene_info[gene1]['start_pos']
    pos2 = gene_info[gene2]['start_pos']
    
    # Different chromosomes = Dispersed or WGD
    if chr1 != chr2:
        if WGD_PEAK1_MIN <= ks <= WGD_PEAK1_MAX or WGD_PEAK2_MIN <= ks <= WGD_PEAK2_MAX:
            return 'WGD', np.nan, abs(pos1 - pos2)
        else:
            return 'Dispersed', np.nan, abs(pos1 - pos2)
    
    # Same chromosome - calculate distances
    bp_distance = abs(pos1 - pos2)
    
    # Calculate gene distance
    if gene1 in gene_order and gene2 in gene_order:
        if gene_order[gene1]['chromosome'] == gene_order[gene2]['chromosome']:
            gene_distance = abs(gene_order[gene1]['order'] - gene_order[gene2]['order'])
        else:
            gene_distance = np.nan
    else:
        gene_distance = np.nan
    
    # Classify
    # TAGs: within 10 genes OR within 100kb
    if (not pd.isna(gene_distance) and gene_distance <= TAG_GENE_DISTANCE) or bp_distance <= TAG_BP_DISTANCE:
        return 'TAG', gene_distance, bp_distance
    
    # Proximal: same chr, 100kb-500kb apart
    elif PROXIMAL_BP_MIN < bp_distance <= PROXIMAL_BP_MAX:
        return 'Proximal', gene_distance, bp_distance
    
    # Dispersed or WGD: same chr, >500kb
    else:
        if WGD_PEAK1_MIN <= ks <= WGD_PEAK1_MAX or WGD_PEAK2_MIN <= ks <= WGD_PEAK2_MAX:
            return 'WGD', gene_distance, bp_distance
        else:
            return 'Dispersed', gene_distance, bp_distance

def classify_all_pairs(ks_df, gene_info, gene_order):
    """Classify all gene pairs."""
    print("\nClassifying gene pairs...")
    
    results = ks_df.apply(lambda row: classify_duplication_type(row, gene_info, gene_order), axis=1)
    
    ks_df['duplication_type'] = [r[0] for r in results]
    ks_df['gene_distance'] = [r[1] for r in results]
    ks_df['bp_distance'] = [r[2] for r in results]
    
    # Count by type
    type_counts = ks_df['duplication_type'].value_counts()
    print("\n  Classification results:")
    for dtype, count in type_counts.items():
        pct = count / len(ks_df) * 100
        print(f"    {dtype}: {count:,} ({pct:.2f}%)")
    
    return ks_df

def create_summary_statistics(ks_df):
    """Create summary statistics by duplication type."""
    print("\nGenerating summary statistics...")
    
    summary_data = []
    
    for dtype in ks_df['duplication_type'].unique():
        subset = ks_df[ks_df['duplication_type'] == dtype]
        
        row = {
            'Duplication_Type': dtype,
            'N_pairs': len(subset),
            'Percentage': f"{len(subset)/len(ks_df)*100:.2f}%",
            'Mean_Ks': subset['ks'].mean(),
            'Median_Ks': subset['ks'].median(),
            'SD_Ks': subset['ks'].std(),
            'Mean_Ka': subset['ka'].mean() if 'ka' in subset.columns else np.nan,
            'Median_Ka_Ks': subset['ka_ks'].median() if 'ka_ks' in subset.columns else np.nan,
        }
        
        # Add distance info for same-chromosome types
        if dtype in ['TAG', 'Proximal']:
            row['Mean_bp_distance'] = subset['bp_distance'].mean()
            row['Median_bp_distance'] = subset['bp_distance'].median()
            if dtype == 'TAG':
                row['Mean_gene_distance'] = subset['gene_distance'].mean()
                row['Median_gene_distance'] = subset['gene_distance'].median()
        
        summary_data.append(row)
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('N_pairs', ascending=False)
    
    # Save to file
    summary_df.to_csv(f"{OUTPUT_DIR}/duplication_type_summary.tsv", sep='\t', index=False, float_format='%.4f')
    print(f"  ✓ Saved summary statistics")
    
    return summary_df

def plot_duplication_types(ks_df, summary_df):
    """Create visualizations of duplication types."""
    print("\nGenerating visualizations...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    fig.suptitle('Gene Duplication Type Classification', fontsize=16, fontweight='bold')
    
    # 1. Pie chart of duplication types
    ax1 = fig.add_subplot(gs[0, 0])
    type_counts = ks_df['duplication_type'].value_counts()
    colors = plt.cm.Set3(range(len(type_counts)))
    wedges, texts, autotexts = ax1.pie(type_counts.values, labels=type_counts.index,
                                         autopct='%1.1f%%', colors=colors, startangle=90)
    for autotext in autotexts:
        autotext.set_fontweight('bold')
        autotext.set_fontsize(9)
    ax1.set_title('Distribution of Duplication Types', fontweight='bold')
    
    # 2. Bar chart of counts
    ax2 = fig.add_subplot(gs[0, 1:])
    type_counts.plot(kind='barh', ax=ax2, color=colors)
    ax2.set_xlabel('Number of gene pairs')
    ax2.set_title('Gene Pair Count by Duplication Type', fontweight='bold')
    ax2.grid(alpha=0.3, axis='x')
    for i, v in enumerate(type_counts.values):
        ax2.text(v + max(type_counts)*0.01, i, f'{v:,}', va='center', fontweight='bold')
    
    # 3. Ks distribution by type (overlaid)
    ax3 = fig.add_subplot(gs[1, :])
    for dtype, color in zip(type_counts.index, colors):
        subset = ks_df[ks_df['duplication_type'] == dtype]['ks']
        ax3.hist(subset, bins=50, alpha=0.5, label=f"{dtype} (n={len(subset)})",
                color=color, density=True)
    ax3.set_xlabel('Ks')
    ax3.set_ylabel('Density')
    ax3.set_title('Ks Distribution by Duplication Type', fontweight='bold')
    ax3.legend()
    ax3.set_xlim(0, 2)
    ax3.grid(alpha=0.3, axis='y')
    
    # Add WGD peak ranges
    ax3.axvspan(WGD_PEAK1_MIN, WGD_PEAK1_MAX, alpha=0.1, color='red', label='WGD Peak 1 (~13 MYA)')
    ax3.axvspan(WGD_PEAK2_MIN, WGD_PEAK2_MAX, alpha=0.1, color='blue', label='WGD Peak 2 (~59 MYA)')
    
    # 4. Box plot of Ks by type
    ax4 = fig.add_subplot(gs[2, 0])
    ks_by_type = [ks_df[ks_df['duplication_type'] == dt]['ks'].values 
                  for dt in type_counts.index]
    bp = ax4.boxplot(ks_by_type, labels=type_counts.index, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    ax4.set_ylabel('Ks')
    ax4.set_title('Ks Distribution (Box Plot)', fontweight='bold')
    ax4.grid(alpha=0.3, axis='y')
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # 5. Ka/Ks by type
    if 'ka_ks' in ks_df.columns:
        ax5 = fig.add_subplot(gs[2, 1])
        kaks_by_type = [ks_df[ks_df['duplication_type'] == dt]['ka_ks'].dropna().values 
                        for dt in type_counts.index]
        bp2 = ax5.boxplot(kaks_by_type, labels=type_counts.index, patch_artist=True)
        for patch, color in zip(bp2['boxes'], colors):
            patch.set_facecolor(color)
        ax5.set_ylabel('Ka/Ks')
        ax5.set_title('Ka/Ks Ratio by Type', fontweight='bold')
        ax5.axhline(1.0, color='red', linestyle='--', alpha=0.5, label='Ka/Ks=1')
        ax5.grid(alpha=0.3, axis='y')
        ax5.legend()
        plt.setp(ax5.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # 6. Summary statistics table
    ax6 = fig.add_subplot(gs[2, 2])
    ax6.axis('tight')
    ax6.axis('off')
    
    table_data = summary_df[['Duplication_Type', 'N_pairs', 'Mean_Ks', 'Median_Ks']].copy()
    table_data['Mean_Ks'] = table_data['Mean_Ks'].apply(lambda x: f'{x:.3f}')
    table_data['Median_Ks'] = table_data['Median_Ks'].apply(lambda x: f'{x:.3f}')
    table_data['N_pairs'] = table_data['N_pairs'].apply(lambda x: f'{x:,}')
    
    table = ax6.table(cellText=table_data.values, colLabels=table_data.columns,
                     cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)
    
    # Style header
    for i in range(len(table_data.columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax6.set_title('Summary Statistics', fontweight='bold', pad=20)
    
    plt.savefig(f"{OUTPUT_DIR}/duplication_types_summary.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{OUTPUT_DIR}/duplication_types_summary.pdf", bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved duplication type visualizations")

def main():
    print("="*60)
    print("DUPLICATION TYPE CLASSIFICATION")
    print("="*60)
    
    # Load data
    gene_info, ks_df = load_data()
    
    # Calculate gene order
    protein_df = pd.read_csv(PROTEIN_INFO)
    gene_order = calculate_gene_order(protein_df)
    
    # Classify pairs
    ks_classified = classify_all_pairs(ks_df, gene_info, gene_order)
    
    # Save classified pairs
    output_file = f"{OUTPUT_DIR}/gene_pairs_classified.tsv"
    ks_classified.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
    print(f"\n  ✓ Saved classified gene pairs to: {output_file}")
    
    # Generate summary statistics
    summary_df = create_summary_statistics(ks_classified)
    
    # Create visualizations
    plot_duplication_types(ks_classified, summary_df)
    
    print("\n" + "="*60)
    print("✓ CLASSIFICATION COMPLETE")
    print("="*60)
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("="*60)

if __name__ == "__main__":
    main()
