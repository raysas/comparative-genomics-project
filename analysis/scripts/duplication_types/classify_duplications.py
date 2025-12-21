#!/usr/bin/env python3
"""
Classify gene duplications into three types with MCScanX validation:
- TAG (Tandem Arrayed Genes): From MCScanX tandem file or detected via proximity
- Proximal: Same chromosome, nearby but not tandem (>10 genes or >100kb, <3Mb)
- WGD (Whole Genome Duplication): Validated against MCScanX collinearity/anchors or Ks peak ranges

Inputs:
  --ks               Ks results TSV (gene1, gene2, ks, ka, ka_ks columns)
  --protein-info     Protein/gene metadata CSV (peptide_id, chromosome, start_pos, end_pos)
  --mcscanx-tandem   MCScanX tandem file (optional, for TAG validation)
  --mcscanx-anchors  MCScanX anchors or collinearity file (optional, for WGD validation)
  --outdir           Output directory (default: analysis/duplication_types)
  
Outputs:
  - gene_pairs_classified.tsv: All pairs with duplication type
  - duplication_summary_stats.tsv: Statistics by type (counts, Ks, Ka/Ks)
  - duplication_ks_distributions.png: Ks histograms by type
  - duplication_ka_ks_boxplots.png: Ka/Ks ratio comparison
  - duplication_type_counts.png: Bar chart of duplication frequencies
  - duplication_chromosome_patterns.png: Structural patterns (intra vs inter-chromosomal)
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150

# Classification parameters
TAG_GENE_DISTANCE = 10  # Maximum genes apart for TAGs
TAG_BP_DISTANCE = 100000  # 100kb
PROXIMAL_BP_MAX = 3000000  # 3Mb
WGD_KS_PEAK1_MIN = 0.05  # Recent WGD (~13 MYA)
WGD_KS_PEAK1_MAX = 0.25
WGD_KS_PEAK2_MIN = 0.40  # Ancient WGD (~59 MYA)
WGD_KS_PEAK2_MAX = 0.65



def load_mcscanx_tandem(tandem_path: Path) -> set:
    """Load MCScanX tandem pairs and return set of (gene1, gene2) tuples."""
    tandem_pairs = set()
    if tandem_path and tandem_path.exists():
        with open(tandem_path) as f:
            for line in f:
                genes = line.strip().split(',')
                if len(genes) >= 2:
                    g1, g2 = genes[0].strip(), genes[1].strip()
                    # Store both directions for lookup
                    tandem_pairs.add((g1, g2))
                    tandem_pairs.add((g2, g1))
        print(f"  ✓ Loaded {len(tandem_pairs)//2} tandem pairs from MCScanX")
    return tandem_pairs


def load_mcscanx_anchors(anchors_path: Path) -> set:
    """Load MCScanX collinearity/anchors and return set of (gene1, gene2) tuples."""
    anchor_pairs = set()
    if anchors_path and anchors_path.exists():
        try:
            # Try as TSV with headers (anchors_with_ks.tsv or glycine.collinearity.anchors.tsv format)
            df = pd.read_csv(anchors_path, sep='\t')
            # Look for gene columns (could be gene1/gene2 or other names)
            gene_cols = [col for col in df.columns if 'gene' in col.lower()]
            if len(gene_cols) >= 2:
                g1_col, g2_col = gene_cols[0], gene_cols[1]
                for _, row in df.iterrows():
                    g1, g2 = str(row[g1_col]), str(row[g2_col])
                    anchor_pairs.add((g1, g2))
                    anchor_pairs.add((g2, g1))
                print(f"  ✓ Loaded {len(anchor_pairs)//2:,} anchor pairs from MCScanX (TSV format)")
                return anchor_pairs
        except Exception as e:
            print(f"  Warning: Failed to parse as TSV: {e}")
        
        # Fallback: try as raw collinearity format (skip headers, parse gene pairs)
        try:
            with open(anchors_path) as f:
                for line in f:
                    if line.startswith('#') or line.startswith('##') or not line.strip():
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        # Collinearity format: ID gene1 gene2 score ...
                        # or block-aligned format:  0-  0:	gene1	gene2	score
                        if ':' in parts[0]:
                            # block-aligned format
                            if len(parts) >= 3:
                                g1, g2 = parts[1], parts[2]
                                anchor_pairs.add((g1, g2))
                                anchor_pairs.add((g2, g1))
                        else:
                            # Simple format: ID gene1 gene2
                            g1, g2 = parts[1], parts[2]
                            anchor_pairs.add((g1, g2))
                            anchor_pairs.add((g2, g1))
            if anchor_pairs:
                print(f"  ✓ Loaded {len(anchor_pairs)//2:,} anchor pairs from MCScanX (raw collinearity format)")
        except Exception as e:
            print(f"  Warning: Failed to parse collinearity file: {e}")
    
    return anchor_pairs


def load_data(ks_path: Path, protein_info_path: Path):
    """Load Ks results and gene/protein metadata."""
    print("Loading data...")
    
    # Load Ks results
    ks_df = pd.read_csv(ks_path, sep='\t')
    print(f"  ✓ Loaded {len(ks_df)} gene pairs with Ks")
    
    # Load protein/gene information
    protein_df = pd.read_csv(protein_info_path)
    print(f"  ✓ Loaded {len(protein_df)} genes")
    
    # Create gene info mapping
    gene_info = protein_df.set_index('peptide_id')[['chromosome', 'start_pos', 'end_pos']].to_dict('index')
    
    return ks_df, protein_df, gene_info


def calculate_gene_order(protein_df):
    """Calculate gene order on each chromosome for proximity detection."""
    print("Calculating gene order on chromosomes...")
    
    gene_order = {}
    for chrom in protein_df['chromosome'].unique():
        chr_genes = protein_df[protein_df['chromosome'] == chrom].sort_values('start_pos')
        for idx, gene_id in enumerate(chr_genes['peptide_id']):
            gene_order[gene_id] = {
                'chromosome': chrom,
                'order': idx,
            }
    
    print(f"  ✓ Calculated order for {len(gene_order)} genes")
    return gene_order


def classify_pair(row, gene_info, gene_order, tandem_pairs, anchor_pairs):
    """Classify a single gene pair into TAG, Proximal, or WGD."""
    gene1 = row['gene1']
    gene2 = row['gene2']
    ks = row['ks']
    
    # Get gene info
    if gene1 not in gene_info or gene2 not in gene_info:
        return 'Unknown', np.nan, np.nan, False, False
    
    chr1 = gene_info[gene1]['chromosome']
    chr2 = gene_info[gene2]['chromosome']
    pos1 = gene_info[gene1]['start_pos']
    pos2 = gene_info[gene2]['start_pos']
    
    bp_distance = abs(pos1 - pos2) if chr1 == chr2 else np.nan
    gene_distance = np.nan
    
    # Calculate gene distance for same chromosome
    if chr1 == chr2 and gene1 in gene_order and gene2 in gene_order:
        gene_distance = abs(gene_order[gene1]['order'] - gene_order[gene2]['order'])
    
    # Check MCScanX validation
    is_mcscanx_tandem = (gene1, gene2) in tandem_pairs or (gene2, gene1) in tandem_pairs
    is_mcscanx_anchor = (gene1, gene2) in anchor_pairs or (gene2, gene1) in anchor_pairs
    
    # Classification logic
    # 1. TAG: ALL proximity-based tandems (MCScanX tandem OR same chr + ≤10 genes apart OR ≤100kb)
    #    This includes both MCScanX-validated and distance-based tandems
    if is_mcscanx_tandem:
        return 'TAG', gene_distance, bp_distance, True, is_mcscanx_anchor
    
    if chr1 == chr2:
        if (not pd.isna(gene_distance) and gene_distance <= TAG_GENE_DISTANCE) or \
           (not pd.isna(bp_distance) and bp_distance <= TAG_BP_DISTANCE):
            return 'TAG', gene_distance, bp_distance, False, is_mcscanx_anchor
        
        # 2. Proximal: Same chr, beyond TAG threshold but < 3Mb
        if not pd.isna(bp_distance) and bp_distance <= PROXIMAL_BP_MAX:
            return 'Proximal', gene_distance, bp_distance, False, is_mcscanx_anchor
    
    # 3. WGD: ONLY MCScanX anchor-validated pairs (filtered, high-confidence WGD)
    #    Do NOT use Ks peaks alone - require synteny evidence
    if is_mcscanx_anchor:
        return 'WGD', gene_distance, bp_distance, is_mcscanx_tandem, True
    
    # 4. Default: Other pairs are Proximal if same chr, otherwise excluded/other
    #    (Inter-chromosomal pairs without MCScanX validation are not classified as WGD)
    if chr1 == chr2:
        return 'Proximal', gene_distance, bp_distance, False, False
    else:
        return 'Other', gene_distance, bp_distance, False, False


def classify_all_pairs(ks_df, gene_info, gene_order, tandem_pairs, anchor_pairs):
    """Classify all gene pairs into duplication types."""
    print("\nClassifying gene pairs...")
    
    results = ks_df.apply(
        lambda row: classify_pair(row, gene_info, gene_order, tandem_pairs, anchor_pairs), 
        axis=1
    )
    
    ks_df['duplication_type'] = [r[0] for r in results]
    ks_df['gene_distance'] = [r[1] for r in results]
    ks_df['bp_distance'] = [r[2] for r in results]
    ks_df['mcscanx_tandem'] = [r[3] for r in results]
    ks_df['mcscanx_anchor'] = [r[4] for r in results]
    
    # Count by type
    type_counts = ks_df['duplication_type'].value_counts()
    print("\n  Classification results:")
    for dtype in ['TAG', 'Proximal', 'WGD', 'Other', 'Unknown']:
        count = type_counts.get(dtype, 0)
        pct = count / len(ks_df) * 100 if len(ks_df) > 0 else 0
        print(f"    {dtype}: {count:,} ({pct:.2f}%)")
    
    # MCScanX validation stats
    if tandem_pairs:
        mcscanx_tag_count = ks_df['mcscanx_tandem'].sum()
        tag_total = (ks_df['duplication_type'] == 'TAG').sum()
        print(f"\n  TAG duplications: {tag_total:,} total (MCScanX validated: {mcscanx_tag_count}, distance-based: {tag_total - mcscanx_tag_count})")
    if anchor_pairs:
        mcscanx_wgd_count = ks_df['mcscanx_anchor'].sum()
        wgd_total = (ks_df['duplication_type'] == 'WGD').sum()
        print(f"  WGD duplications: {wgd_total:,} total (all MCScanX anchor-validated)")
    
    return ks_df




def create_summary_statistics(ks_df, outdir: Path):
    """Generate summary statistics by duplication type."""
    print("\nGenerating summary statistics...")
    
    summary_rows = []
    
    for dtype in ['TAG', 'Proximal', 'WGD']:
        subset = ks_df[ks_df['duplication_type'] == dtype]
        if len(subset) == 0:
            continue
        
        row = {
            'Duplication_Type': dtype,
            'N_pairs': len(subset),
            'Percentage': f"{len(subset)/len(ks_df)*100:.2f}%",
            'Mean_Ks': subset['ks'].mean(),
            'Median_Ks': subset['ks'].median(),
            'SD_Ks': subset['ks'].std(),
            'Min_Ks': subset['ks'].min(),
            'Max_Ks': subset['ks'].max(),
        }
        
        # Ka/Ks if available
        if 'ka_ks' in subset.columns:
            ka_ks_valid = subset['ka_ks'].dropna()
            if len(ka_ks_valid) > 0:
                row['Mean_Ka_Ks'] = ka_ks_valid.mean()
                row['Median_Ka_Ks'] = ka_ks_valid.median()
        
        # Add distance info for TAG/Proximal
        if dtype in ['TAG', 'Proximal']:
            bp_dist_valid = subset['bp_distance'].dropna()
            if len(bp_dist_valid) > 0:
                row['Mean_bp_distance'] = bp_dist_valid.mean()
                row['Median_bp_distance'] = bp_dist_valid.median()
            
            if dtype == 'TAG':
                gene_dist_valid = subset['gene_distance'].dropna()
                if len(gene_dist_valid) > 0:
                    row['Mean_gene_distance'] = gene_dist_valid.mean()
                    row['Median_gene_distance'] = gene_dist_valid.median()
        
        # MCScanX validation counts
        if 'mcscanx_tandem' in subset.columns:
            row['MCScanX_tandem_validated'] = subset['mcscanx_tandem'].sum()
        if 'mcscanx_anchor' in subset.columns:
            row['MCScanX_anchor_validated'] = subset['mcscanx_anchor'].sum()
        
        summary_rows.append(row)
    
    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values('N_pairs', ascending=False)
    
    # Save
    summary_path = outdir / 'duplication_summary_stats.tsv'
    summary_df.to_csv(summary_path, sep='\t', index=False, float_format='%.4f')
    print(f"  ✓ Saved summary statistics to {summary_path}")
    
    return summary_df


def plot_ks_distributions(ks_df, outdir: Path):
    """Plot Ks distributions by duplication type (TAG, Proximal, WGD only)."""
    print("Generating Ks distribution plots...")
    
    # Filter to only the three main duplication types for analysis
    ks_df_filtered = ks_df[ks_df['duplication_type'].isin(['TAG', 'Proximal', 'WGD'])].copy()
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Ks Distributions by Duplication Type', fontsize=14, fontweight='bold')
    
    colors = {'TAG': '#e74c3c', 'Proximal': '#3498db', 'WGD': '#2ecc71'}
    
    # 1. Overlaid histograms
    ax = axes[0, 0]
    for dtype in ['TAG', 'Proximal', 'WGD']:
        subset = ks_df_filtered[ks_df_filtered['duplication_type'] == dtype]['ks']
        if len(subset) > 0:
            ax.hist(subset, bins=60, alpha=0.6, label=f'{dtype} (n={len(subset):,})',
                   color=colors.get(dtype, 'gray'), density=True, range=(0, 2))
    ax.set_xlabel('Ks')
    ax.set_ylabel('Density')
    ax.set_title('Overlaid Ks Distributions')
    ax.legend(loc='upper right')
    ax.grid(alpha=0.3, axis='y')
    ax.set_xlim(0, 2)
    
    # Add WGD peak ranges
    ax.axvspan(WGD_KS_PEAK1_MIN, WGD_KS_PEAK1_MAX, alpha=0.15, color='orange', 
              label='WGD Peak 1', zorder=0)
    ax.axvspan(WGD_KS_PEAK2_MIN, WGD_KS_PEAK2_MAX, alpha=0.15, color='purple', 
              label='WGD Peak 2', zorder=0)
    
    # 2. Separate subplots
    ax = axes[0, 1]
    for i, dtype in enumerate(['TAG', 'Proximal', 'WGD']):
        subset = ks_df_filtered[ks_df_filtered['duplication_type'] == dtype]['ks']
        if len(subset) > 0:
            ax.hist(subset, bins=50, alpha=0.7, label=dtype,
                   color=colors.get(dtype, 'gray'), histtype='step', linewidth=2, range=(0, 2))
    ax.set_xlabel('Ks')
    ax.set_ylabel('Frequency')
    ax.set_title('Ks Distributions (line plot)')
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 2)
    
    # 3. Box plot
    ax = axes[1, 0]
    ks_by_type = []
    labels = []
    box_colors = []
    for dtype in ['TAG', 'Proximal', 'WGD']:
        subset = ks_df_filtered[ks_df_filtered['duplication_type'] == dtype]['ks'].dropna()
        if len(subset) > 0:
            ks_by_type.append(subset.values)
            labels.append(f'{dtype}\n(n={len(subset)})')
            box_colors.append(colors.get(dtype, 'gray'))
    
    bp = ax.boxplot(ks_by_type, labels=labels, patch_artist=True, widths=0.6)
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel('Ks')
    ax.set_title('Ks Distributions (box plot)')
    ax.grid(alpha=0.3, axis='y')
    
    # 4. Violin plot
    ax = axes[1, 1]
    plot_data = []
    for dtype in ['TAG', 'Proximal', 'WGD']:
        subset = ks_df_filtered[ks_df_filtered['duplication_type'] == dtype][['ks']].copy()
        subset['Type'] = dtype
        plot_data.append(subset)
    if plot_data:
        plot_df = pd.concat(plot_data)
        sns.violinplot(data=plot_df, x='Type', y='ks', palette=colors, ax=ax, inner='box')
        ax.set_ylabel('Ks')
        ax.set_xlabel('')
        ax.set_title('Ks Distributions (violin plot)')
        ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    outpath = outdir / 'duplication_ks_distributions.png'
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'duplication_ks_distributions.pdf', bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved Ks distributions to {outpath}")


def plot_ka_ks_ratios(ks_df, outdir: Path):
    """Plot Ka/Ks ratios by duplication type (TAG, Proximal, WGD only)."""
    print("Generating Ka/Ks ratio plots...")
    
    if 'ka_ks' not in ks_df.columns:
        print("  ⚠ Ka/Ks column not found, skipping Ka/Ks plots")
        return
    
    # Filter to only the three main duplication types
    ks_df_filtered = ks_df[ks_df['duplication_type'].isin(['TAG', 'Proximal', 'WGD'])].copy()
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Ka/Ks Ratios by Duplication Type', fontsize=14, fontweight='bold')
    
    colors = {'TAG': '#e74c3c', 'Proximal': '#3498db', 'WGD': '#2ecc71'}
    
    # 1. Box plot
    ax = axes[0]
    kaks_by_type = []
    labels = []
    box_colors = []
    for dtype in ['TAG', 'Proximal', 'WGD']:
        subset = ks_df_filtered[ks_df_filtered['duplication_type'] == dtype]['ka_ks'].dropna()
        # Filter extreme outliers
        subset = subset[(subset >= 0) & (subset <= 3)]
        if len(subset) > 0:
            kaks_by_type.append(subset.values)
            labels.append(f'{dtype}\n(n={len(subset)})')
            box_colors.append(colors.get(dtype, 'gray'))
    
    if kaks_by_type:
        bp = ax.boxplot(kaks_by_type, labels=labels, patch_artist=True, widths=0.6)
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.axhline(1.0, color='red', linestyle='--', alpha=0.6, linewidth=1.5, label='Ka/Ks=1 (neutral)')
        ax.set_ylabel('Ka/Ks')
        ax.set_title('Ka/Ks Ratio Distribution')
        ax.grid(alpha=0.3, axis='y')
        ax.legend()
        ax.set_ylim(0, 2.5)
    
    # 2. Violin plot
    ax = axes[1]
    plot_data = []
    for dtype in ['TAG', 'Proximal', 'WGD']:
        subset = ks_df_filtered[ks_df_filtered['duplication_type'] == dtype][['ka_ks']].copy()
        subset = subset[(subset['ka_ks'] >= 0) & (subset['ka_ks'] <= 3)].dropna()
        subset['Type'] = dtype
        plot_data.append(subset)
    
    if plot_data:
        plot_df = pd.concat(plot_data)
        if len(plot_df) > 0:
            sns.violinplot(data=plot_df, x='Type', y='ka_ks', palette=colors, ax=ax, inner='quartile')
            ax.axhline(1.0, color='red', linestyle='--', alpha=0.6, linewidth=1.5, label='Ka/Ks=1')
            ax.set_ylabel('Ka/Ks')
            ax.set_xlabel('')
            ax.set_title('Ka/Ks Ratio (violin plot)')
            ax.grid(alpha=0.3, axis='y')
            ax.legend()
            ax.set_ylim(0, 2.5)
    
    plt.tight_layout()
    outpath = outdir / 'duplication_ka_ks_boxplots.png'
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'duplication_ka_ks_boxplots.pdf', bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved Ka/Ks plots to {outpath}")


def plot_type_counts(ks_df, outdir: Path):
    """Plot duplication type frequencies (TAG, Proximal, WGD only)."""
    print("Generating duplication type count plots...")
    
    # Filter to only the three main duplication types
    ks_df_filtered = ks_df[ks_df['duplication_type'].isin(['TAG', 'Proximal', 'WGD'])].copy()
    total_analyzed = len(ks_df_filtered)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Duplication Type Frequencies', fontsize=14, fontweight='bold')
    
    colors = {'TAG': '#e74c3c', 'Proximal': '#3498db', 'WGD': '#2ecc71'}
    type_counts = ks_df_filtered['duplication_type'].value_counts()
    
    # Reorder to TAG, Proximal, WGD
    ordered_types = [t for t in ['TAG', 'Proximal', 'WGD'] if t in type_counts.index]
    type_counts = type_counts.loc[ordered_types]
    
    # 1. Bar chart
    ax = axes[0]
    bar_colors = [colors.get(t, 'gray') for t in type_counts.index]
    bars = ax.bar(range(len(type_counts)), type_counts.values, color=bar_colors, alpha=0.8, edgecolor='black')
    ax.set_xticks(range(len(type_counts)))
    ax.set_xticklabels(type_counts.index)
    ax.set_ylabel('Number of gene pairs')
    ax.set_title('Duplication Type Counts')
    ax.grid(alpha=0.3, axis='y')
    
    # Add count labels on bars
    for i, (count, bar) in enumerate(zip(type_counts.values, bars)):
        pct = count / total_analyzed * 100
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(type_counts)*0.01,
               f'{count:,}\n({pct:.1f}%)', ha='center', va='bottom', fontweight='bold', fontsize=10)
    
    # 2. Pie chart
    ax = axes[1]
    wedges, texts, autotexts = ax.pie(type_counts.values, labels=type_counts.index,
                                        autopct='%1.1f%%', colors=bar_colors, startangle=90,
                                        textprops={'fontsize': 11, 'fontweight': 'bold'})
    ax.set_title('Duplication Type Proportions')
    
    plt.tight_layout()
    outpath = outdir / 'duplication_type_counts.png'
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'duplication_type_counts.pdf', bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved type count plots to {outpath}")


def plot_chromosome_patterns(ks_df, gene_info, outdir: Path):
    """Analyze and plot structural/chromosomal patterns (TAG, Proximal, WGD only)."""
    print("Generating chromosome pattern analysis...")
    
    # Filter to only the three main duplication types
    ks_df_filtered = ks_df[ks_df['duplication_type'].isin(['TAG', 'Proximal', 'WGD'])].copy()
    
    # Add chromosome info to pairs
    ks_df_chr = ks_df_filtered.copy()
    ks_df_chr['chr1'] = ks_df_chr['gene1'].map(lambda g: gene_info.get(g, {}).get('chromosome', None))
    ks_df_chr['chr2'] = ks_df_chr['gene2'].map(lambda g: gene_info.get(g, {}).get('chromosome', None))
    ks_df_chr['same_chr'] = ks_df_chr['chr1'] == ks_df_chr['chr2']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Chromosomal Distribution Patterns', fontsize=14, fontweight='bold')
    
    colors = {'TAG': '#e74c3c', 'Proximal': '#3498db', 'WGD': '#2ecc71'}
    
    # 1. Intra vs inter-chromosomal by type
    ax = axes[0, 0]
    pattern_data = []
    for dtype in ['TAG', 'Proximal', 'WGD']:
        subset = ks_df_chr[ks_df_chr['duplication_type'] == dtype]
        if len(subset) > 0:
            intra = subset['same_chr'].sum()
            inter = len(subset) - intra
            pattern_data.append({'Type': dtype, 'Intra-chromosomal': intra, 'Inter-chromosomal': inter})
    
    if pattern_data:
        pattern_df = pd.DataFrame(pattern_data).set_index('Type')
        pattern_df.plot(kind='bar', stacked=False, ax=ax, color=['#95a5a6', '#34495e'], alpha=0.8)
        ax.set_ylabel('Number of pairs')
        ax.set_title('Intra- vs Inter-chromosomal Duplications')
        ax.legend(loc='upper right')
        ax.grid(alpha=0.3, axis='y')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    
    # 2. Percentage stacked bar
    ax = axes[0, 1]
    if pattern_data:
        pattern_df_pct = pattern_df.div(pattern_df.sum(axis=1), axis=0) * 100
        pattern_df_pct.plot(kind='bar', stacked=True, ax=ax, color=['#95a5a6', '#34495e'], alpha=0.8)
        ax.set_ylabel('Percentage')
        ax.set_title('Intra- vs Inter-chromosomal (Percentage)')
        ax.legend(loc='upper right')
        ax.grid(alpha=0.3, axis='y')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
    
    # 3. Chromosome pair frequency heatmap for WGD
    ax = axes[1, 0]
    wgd_subset = ks_df_chr[ks_df_chr['duplication_type'] == 'WGD']
    if len(wgd_subset) > 0:
        chr_pair_counts = defaultdict(int)
        for _, row in wgd_subset.iterrows():
            c1, c2 = row['chr1'], row['chr2']
            if pd.notna(c1) and pd.notna(c2):
                pair = tuple(sorted([str(c1), str(c2)]))
                chr_pair_counts[pair] += 1
        
        # Top chromosome pairs
        top_pairs = sorted(chr_pair_counts.items(), key=lambda x: x[1], reverse=True)[:15]
        if top_pairs:
            pair_labels = [f"{p[0]}-{p[1]}" for p, _ in top_pairs]
            pair_counts = [c for _, c in top_pairs]
            ax.barh(range(len(pair_labels)), pair_counts, color=colors['WGD'], alpha=0.8)
            ax.set_yticks(range(len(pair_labels)))
            ax.set_yticklabels(pair_labels, fontsize=9)
            ax.set_xlabel('Number of WGD pairs')
            ax.set_title('Top 15 Chromosome Pairs (WGD)')
            ax.grid(alpha=0.3, axis='x')
    
    # 4. Gene distance distribution for TAGs
    ax = axes[1, 1]
    tag_subset = ks_df_chr[ks_df_chr['duplication_type'] == 'TAG']
    gene_dist = tag_subset['gene_distance'].dropna()
    if len(gene_dist) > 0:
        ax.hist(gene_dist, bins=range(0, min(int(gene_dist.max())+2, 21)), 
               color=colors['TAG'], alpha=0.7, edgecolor='black')
        ax.set_xlabel('Genes apart')
        ax.set_ylabel('Frequency')
        ax.set_title(f'TAG Gene Distance Distribution (n={len(gene_dist)})')
        ax.grid(alpha=0.3, axis='y')
        ax.set_xlim(0, 20)
    
    plt.tight_layout()
    outpath = outdir / 'duplication_chromosome_patterns.png'
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'duplication_chromosome_patterns.pdf', bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved chromosome pattern plots to {outpath}")




def main():
    parser = argparse.ArgumentParser(
        description='Classify gene duplications into TAG, Proximal, and WGD types'
    )
    parser.add_argument('--ks', type=Path, default=Path('pipeline_2/ks_results_filtered.tsv'),
                       help='Ks results TSV file (gene1, gene2, ks, ka, ka_ks)')
    parser.add_argument('--protein-info', type=Path, default=Path('data/protein_info_longest.csv'),
                       help='Protein/gene metadata CSV')
    parser.add_argument('--mcscanx-tandem', type=Path, default=Path('output/mcscanx/glycine.tandem'),
                       help='MCScanX tandem file (for TAG validation)')
    parser.add_argument('--mcscanx-anchors', type=Path, default=Path('output/mcscanx/glycine.collinearity.anchors.tsv'),
                       help='MCScanX anchors/collinearity file (for WGD validation)')
    parser.add_argument('--outdir', type=Path, default=Path('analysis/duplication_types'),
                       help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    args.outdir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("GENE DUPLICATION CLASSIFICATION: TAG / Proximal / WGD")
    print("="*70)
    
    # Load MCScanX data if provided
    tandem_pairs = set()
    anchor_pairs = set()
    
    if args.mcscanx_tandem:
        tandem_pairs = load_mcscanx_tandem(args.mcscanx_tandem)
    
    if args.mcscanx_anchors:
        anchor_pairs = load_mcscanx_anchors(args.mcscanx_anchors)
    
    # Load main data
    ks_df, protein_df, gene_info = load_data(args.ks, args.protein_info)
    
    # Calculate gene order
    gene_order = calculate_gene_order(protein_df)
    
    # Classify all pairs
    ks_classified = classify_all_pairs(ks_df, gene_info, gene_order, tandem_pairs, anchor_pairs)
    
    # Save classified pairs
    output_file = args.outdir / 'gene_pairs_classified.tsv'
    ks_classified.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
    print(f"\n  ✓ Saved classified gene pairs to: {output_file}")
    
    # Generate summary statistics
    summary_df = create_summary_statistics(ks_classified, args.outdir)
    
    # Generate all plots
    plot_ks_distributions(ks_classified, args.outdir)
    plot_ka_ks_ratios(ks_classified, args.outdir)
    plot_type_counts(ks_classified, args.outdir)
    plot_chromosome_patterns(ks_classified, gene_info, args.outdir)
    
    print("\n" + "="*70)
    print("✓ CLASSIFICATION COMPLETE")
    print("="*70)
    print(f"\nOutputs saved to: {args.outdir}")
    print("\nGenerated files:")
    print("  - gene_pairs_classified.tsv")
    print("  - duplication_summary_stats.tsv")
    print("  - duplication_ks_distributions.png/pdf")
    print("  - duplication_ka_ks_boxplots.png/pdf")
    print("  - duplication_type_counts.png/pdf")
    print("  - duplication_chromosome_patterns.png/pdf")
    print("="*70)


if __name__ == '__main__':
    main()

