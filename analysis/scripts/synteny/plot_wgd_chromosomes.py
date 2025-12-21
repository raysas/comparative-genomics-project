#!/usr/bin/env python3
"""
Plot chromosome–chromosome synteny by WGD and Ks mixture context.

Inputs (from analyze_wgd_chromosomes.py outputs):
  --in-dir  directory containing:
            - chromosome_pairs_wgd.tsv
            - wgd_mixture_stats.tsv
  --out-dir directory to write figures

Outputs:
  - chr_chr_heatmap_by_wgd.png / .pdf
  - ks_mixture_summary.png / .pdf
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("white")


def _numeric_sort_chr(chrs):
    def key(c):
        try:
            return int(str(c).lstrip('chr').lstrip('Chr'))
        except Exception:
            return 9999
    return sorted(chrs, key=key)

def plot_chr_chr_heatmap(pairs_tsv: Path, out_dir: Path, mixture_tsv: Path | None = None, normalize_rows: bool = False):
    df = pd.read_csv(pairs_tsv, sep='\t')
    if df.empty:
        raise RuntimeError("chromosome_pairs_wgd.tsv is empty")
    # Ensure categories sorted numerically where possible
    df['chr1'] = pd.Categorical(df['chr1'], _numeric_sort_chr(df['chr1'].unique()))
    df['chr2'] = pd.Categorical(df['chr2'], _numeric_sort_chr(df['chr2'].unique()))

    # Facet by WGD label
    wgd_labels = ['WGD1_lowKs','WGD2_highKs']
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    cmap = sns.cm.rocket
    # Optional: read mixture stats for title annotation
    mix_annot = {}
    if mixture_tsv is not None and mixture_tsv.exists():
        ms = pd.read_csv(mixture_tsv, sep='\t')
        if set(['component','mean','sd','weight']).issubset(ms.columns):
            for _, r in ms.iterrows():
                mix_annot[str(r['component'])] = (float(r['mean']), float(r['weight']))
    for i, wgd in enumerate(wgd_labels):
        sub = df[df['wgd_label'] == wgd]
        pivot = sub.pivot_table(index='chr1', columns='chr2', values='n_anchors', aggfunc='sum', fill_value=0)
        if normalize_rows and not pivot.empty:
            row_sums = pivot.sum(axis=1).replace(0, np.nan)
            pivot = pivot.div(row_sums, axis=0).fillna(0.0)
        ax = axes[i]
        sns.heatmap(pivot, ax=ax, cmap=cmap, cbar=i==1, annot=False, linewidths=0.2)
        # Top-match markers per row
        if not pivot.empty:
            for r_idx, row in enumerate(pivot.values):
                if normalize_rows:
                    j = int(np.nanargmax(row))
                else:
                    j = int(np.argmax(row))
                ax.scatter(j + 0.5, r_idx + 0.5, s=20, color='#7A0177' if wgd=='WGD2_highKs' else '#084594', zorder=3)
        # Title with optional mixture annotation
        title = f"{wgd} ({'fraction per chr1' if normalize_rows else 'anchor counts'})"
        if wgd in mix_annot:
            mu, wt = mix_annot[wgd]
            title += f"\nμ={mu:.2f}, weight={wt:.2f}"
        ax.set_title(title, fontweight='bold', pad=12)
        ax.set_xlabel('chr2')
        ax.set_ylabel('chr1')
        # rotate labels for readability
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    (out_dir / ('chr_chr_heatmap_by_wgd_normalized.png' if normalize_rows else 'chr_chr_heatmap_by_wgd_counts.png')).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / ('chr_chr_heatmap_by_wgd_normalized.png' if normalize_rows else 'chr_chr_heatmap_by_wgd_counts.png'), dpi=300)
    fig.savefig(out_dir / ('chr_chr_heatmap_by_wgd_normalized.pdf' if normalize_rows else 'chr_chr_heatmap_by_wgd_counts.pdf'))
    plt.close(fig)


def plot_ks_mixture(mixture_tsv: Path, out_dir: Path):
    stats = pd.read_csv(mixture_tsv, sep='\t')
    # Expect columns: component, mean, sd, weight
    if not set(['component','mean','sd','weight']).issubset(stats.columns):
        # Fallback: mixture not used; skip
        return
    x = np.linspace(0, 2, 400)
    # Construct component curves (scaled by weights)
    from scipy.stats import norm
    comp_curves = []
    for _, r in stats.iterrows():
        comp = r['component']
        mu = float(r['mean']); sd = float(r['sd']); w = float(r['weight'])
        y = w * norm.pdf(x, mu, sd)
        comp_curves.append((comp, x, y, mu))
    fig, ax = plt.subplots(figsize=(10, 5))
    # Plot components
    colors = {'WGD1_lowKs':'#084594','WGD2_highKs':'#7A0177'}
    for comp, xs, ys, mu in comp_curves:
        ax.plot(xs, ys, color=colors.get(comp, '#333333'), lw=2, label=f"{comp} (μ={mu:.2f})")
        ax.axvline(mu, color=colors.get(comp, '#333333'), lw=1, ls=':')
    ax.set_xlabel('Ks')
    ax.set_ylabel('Density (scaled)')
    ax.set_title('Gaussian Mixture Components (Ks)', fontweight='bold', pad=12)
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / 'ks_mixture_summary.png', dpi=300)
    fig.savefig(out_dir / 'ks_mixture_summary.pdf')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Plot WGD synteny summaries for slides')
    ap.add_argument('--in-dir', required=True, help='Input dir with chromosome_pairs_wgd.tsv and wgd_mixture_stats.tsv')
    ap.add_argument('--out-dir', required=True, help='Output dir for figures')
    ap.add_argument('--normalize', action='store_true', help='Normalize rows to fractions per chr1')
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_chr_chr_heatmap(in_dir / 'chromosome_pairs_wgd.tsv', out_dir, in_dir / 'wgd_mixture_stats.tsv', normalize_rows=args.normalize)
    # Also render the other variant for convenience
    if not args.normalize:
        plot_chr_chr_heatmap(in_dir / 'chromosome_pairs_wgd.tsv', out_dir, in_dir / 'wgd_mixture_stats.tsv', normalize_rows=True)
    plot_ks_mixture(in_dir / 'wgd_mixture_stats.tsv', out_dir)
    print(f"Figures written to: {out_dir}")


if __name__ == '__main__':
    main()
