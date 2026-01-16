#!/usr/bin/env python3
"""
Clean, presentation-ready summaries of WGD chromosome matches.

Inputs:
  --in    TSV from analyze_wgd_chromosomes.py (chromosome_summary_wgd.tsv)
          columns: chr1,wgd_label,chr2,top_count,top_median_ks
  --out   Output directory for figures

Outputs:
  - matched_bars_per_chr.(png|pdf): side-by-side bars per chr1 for WGD1/WGD2 showing top chr2 counts
  - ranked_pairs_by_wgd.(png|pdf): compact table-style figure listing top matches per WGD
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

WGD_COLORS = {
    'WGD1_lowKs': '#084594',
    'WGD2_highKs': '#7A0177',
}


def numeric_chr_sort(chrs):
    def key(c):
        s = str(c)
        s = s.lstrip('chr').lstrip('Chr')
        try:
            return int(s)
        except Exception:
            return 9999
    return sorted(chrs, key=key)


def plot_matched_bars(summary_path: Path, out_dir: Path):
    df = pd.read_csv(summary_path, sep='\t')
    if df.empty:
        raise RuntimeError('Input summary is empty')
    # Ensure consistent ordering
    chr_order = numeric_chr_sort(df['chr1'].unique())
    wgd_order = ['WGD1_lowKs','WGD2_highKs']

    # Build per-chr1 records for both WGDs
    rows = []
    for chr1 in chr_order:
        for wgd in wgd_order:
            sub = df[(df['chr1'] == chr1) & (df['wgd_label'] == wgd)]
            if not sub.empty:
                r = sub.iloc[0]
                rows.append({'chr1': chr1, 'wgd_label': wgd, 'chr2': r['chr2'], 'count': int(r['top_count'])})
            else:
                rows.append({'chr1': chr1, 'wgd_label': wgd, 'chr2': '', 'count': 0})
    plot_df = pd.DataFrame(rows)

    # Bar plot
    fig, ax = plt.subplots(figsize=(14, 7))
    x = np.arange(len(chr_order))
    width = 0.42
    wgd1 = plot_df[plot_df['wgd_label'] == 'WGD1_lowKs']
    wgd2 = plot_df[plot_df['wgd_label'] == 'WGD2_highKs']
    ax.bar(x - width/2, wgd1['count'].to_numpy(), width, color=WGD_COLORS['WGD1_lowKs'], label='WGD1 (low Ks)')
    ax.bar(x + width/2, wgd2['count'].to_numpy(), width, color=WGD_COLORS['WGD2_highKs'], label='WGD2 (high Ks)')

    # Annotate chr2 labels on bars (small text)
    for i, (c1, lab) in enumerate(zip(wgd1['chr1'], wgd1['chr2'])):
        val = wgd1['count'].iloc[i]
        if val > 0 and lab:
            ax.text(i - width/2, val + max(1, val*0.02), str(lab), ha='center', va='bottom', fontsize=8, color=WGD_COLORS['WGD1_lowKs'])
    for i, (c1, lab) in enumerate(zip(wgd2['chr1'], wgd2['chr2'])):
        val = wgd2['count'].iloc[i]
        if val > 0 and lab:
            ax.text(i + width/2, val + max(1, val*0.02), str(lab), ha='center', va='bottom', fontsize=8, color=WGD_COLORS['WGD2_highKs'])

    ax.set_xticks(x)
    ax.set_xticklabels(chr_order, rotation=0)
    ax.set_ylabel('Top-match anchor count')
    ax.set_xlabel('chr1')
    ax.set_title('Top chromosome matches per chr1 (bars: counts; labels: chr2)', fontweight='bold', pad=14)
    ax.legend(loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / 'matched_bars_per_chr.png', dpi=300)
    fig.savefig(out_dir / 'matched_bars_per_chr.pdf')
    plt.close(fig)


def plot_ranked_pairs(summary_path: Path, out_dir: Path):
    df = pd.read_csv(summary_path, sep='\t')
    if df.empty:
        raise RuntimeError('Input summary is empty')
    # Order by chr1 then WGD
    df['chr1'] = pd.Categorical(df['chr1'], numeric_chr_sort(df['chr1'].unique()))
    df['wgd_label'] = pd.Categorical(df['wgd_label'], ['WGD1_lowKs','WGD2_highKs'])
    df = df.sort_values(['chr1','wgd_label']).reset_index(drop=True)
    # Build a compact table figure using matplotlib
    cols = ['chr1','wgd_label','chr2','top_count','top_median_ks']
    display = df[cols].copy()
    display['wgd_label'] = display['wgd_label'].map({'WGD1_lowKs':'WGD1','WGD2_highKs':'WGD2'})

    fig = plt.figure(figsize=(12, min(12, 0.35*len(display)+2)))
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.axis('off')

    # Build table data
    cell_text = display.values.tolist()
    col_labels = ['chr1','WGD','chr2 (top match)','anchors','median Ks']
    table = ax.table(cellText=cell_text, colLabels=col_labels, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.3)

    # Color WGD column cells
    for i in range(len(cell_text)):
        wgd = display.loc[i, 'wgd_label']
        color = WGD_COLORS['WGD1_lowKs'] if wgd == 'WGD1' else WGD_COLORS['WGD2_highKs']
        # WGD column index = 1
        table[(i+1, 1)].set_facecolor(color+'20')  # translucent
        table[(i+1, 1)].get_text().set_color(color)

    ax.set_title('Ranked top matches per chromosome and WGD', fontweight='bold', pad=14)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / 'ranked_pairs_by_wgd.png', dpi=300)
    fig.savefig(out_dir / 'ranked_pairs_by_wgd.pdf')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Clean WGD match plots for slides')
    ap.add_argument('--in', required=True, dest='infile', help='chromosome_summary_wgd.tsv')
    ap.add_argument('--out', required=True, dest='outdir', help='Output directory for figures')
    args = ap.parse_args()
    in_path = Path(args.infile)
    out_dir = Path(args.outdir)

    plot_matched_bars(in_path, out_dir)
    plot_ranked_pairs(in_path, out_dir)
    print(f"Wrote figures to {out_dir}")


if __name__ == '__main__':
    main()
