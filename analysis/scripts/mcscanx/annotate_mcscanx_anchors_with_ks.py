#!/usr/bin/env python3
"""
Parse MCScanX .collinearity anchors, join with Ks results, and plot.

Inputs:
  --col   output/mcscanx/soybean.collinearity
  --ks    output/ks_results/ks_results_filtered.tsv (columns: gene1, gene2, ks, ...)
  --out   output/mcscanx/soybean_anchors_output/
"""

import argparse
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

ANCHOR_LINE = re.compile(r"^\s*\d+-\s+\d+:\s+(\S+)\s+(\S+)\s+([eE0-9\-\.]+)")

def parse_collinearity(col_path: str) -> pd.DataFrame:
    anchors = []
    block_id = None
    block_score = None
    with open(col_path) as f:
        for line in f:
            if line.startswith("## Alignment"):
                # Header example: ## Alignment 1: score= 34.0 e_value= 1.1e-05 N= 3 ...
                block_id = int(re.search(r"Alignment\s+(\d+)", line).group(1))
                mscore = re.search(r"score\s*=\s*([0-9\.]+)", line)
                block_score = float(mscore.group(1)) if mscore else np.nan
                continue
            m = ANCHOR_LINE.match(line)
            if m:
                g1, g2, evalue = m.group(1), m.group(2), m.group(3)
                try:
                    ev = float(evalue)
                except Exception:
                    # Try to fix malformed scientific notation (e.g., '3e' -> '3e0')
                    if re.match(r'^[0-9]+e$', evalue):
                        ev = float(evalue + '0')
                    else:
                        ev = np.nan
                anchors.append({
                    'block': block_id,
                    'block_score': block_score,
                    'gene1': g1,
                    'gene2': g2,
                    'anchor_evalue': ev
                })
    # save as file
    anchors_df = pd.DataFrame(anchors)
    anchors_df.to_csv(col_path + ".anchors.tsv", sep='\t', index=False)
    return anchors_df


def load_ks(ks_path: str) -> pd.DataFrame:
    ks = pd.read_csv(ks_path, sep='\t')
    # Normalize columns
    ks = ks.rename(columns={c: c.lower() for c in ks.columns})
    return ks


def join_anchors_ks(anchors: pd.DataFrame, ks: pd.DataFrame) -> pd.DataFrame:
    if anchors.empty:
        return anchors
    # Prepare for two-way join
    a = anchors.copy()
    k1 = ks.copy()
    k2 = ks.copy()
    # k1: gene1-gene2 as-is
    j1 = a.merge(k1, left_on=['gene1','gene2'], right_on=['gene1','gene2'], how='left')
    # k2: swap
    j2 = a.merge(k2, left_on=['gene1','gene2'], right_on=['gene2','gene1'], how='left', suffixes=(None,'_sw'))

    # Prefer j1 values; fill from j2 where missing
    def coalesce(s1, s2):
        return s1.where(s1.notna(), s2)

    out = j1.copy()
    for col in ['ks','ka','ka_ks','length','family','status']:
        c2 = col if col in j2.columns else f"{col}_sw"
        if col in out.columns and c2 in j2.columns:
            out[col] = coalesce(out[col], j2[c2])
        elif c2 in j2.columns:
            out[col] = j2[c2]
    return out


def plot_anchors_ks(df: pd.DataFrame, out_dir: str):
    filtered = df[df['ks'].between(0, 2, inclusive='both')]

    plt.figure(figsize=(10,6))
    plt.hist(filtered['ks'].dropna(), bins=100, color='steelblue', alpha=0.8, density=True)
    plt.xlabel('Ks (anchors)')
    plt.ylabel('Density')
    plt.title('Ks Distribution for MCScanX Anchor Pairs')
    plt.tight_layout()
    plt.savefig(f"{out_dir}/ks_hist.png", dpi=300)
    plt.savefig(f"{out_dir}/ks_hist.pdf")
    plt.close()

    # Blocks by mean Ks
    block_stats = filtered.groupby('block')['ks'].agg(['count','median','mean']).reset_index()
    block_stats = block_stats.sort_values('median')

    plt.figure(figsize=(12,6))
    plt.scatter(block_stats['block'], block_stats['median'], s=10, alpha=0.7)
    plt.axhspan(0.1, 0.3, color='red', alpha=0.1, label='WGD ~13MYA')
    plt.axhspan(0.7, 1.2, color='blue', alpha=0.1, label='WGD ~59MYA')
    plt.xlabel('Syntenic block ID')
    plt.ylabel('Median Ks per block')
    plt.title('Median Ks per Syntenic Block (anchors)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_dir}/block_median_ks.png", dpi=300)
    plt.savefig(f"{out_dir}/block_median_ks.pdf")
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--col', required=True, help='MCScanX .collinearity file')
    ap.add_argument('--ks', required=True, help='Ks results TSV (gene1, gene2, ks, ...)')
    ap.add_argument('--out', required=True, help='Output directory for TSV and plots')
    
    args = ap.parse_args()
        
    print("Parsing collinearity file...")
    anchors = parse_collinearity(args.col)
    print(f"  Anchors parsed: {len(anchors)}")

    print("Loading Ks results...")
    ks = load_ks(args.ks)
    print(f"  Ks pairs loaded: {len(ks)}")

    print("Joining anchors with Ks...")
    joined = join_anchors_ks(anchors, ks)
    print(f"  Anchors with Ks annotated: {joined['ks'].notna().sum()} (of {len(joined)})")

    # if output directory does not exist, create it
    Path(args.out).mkdir(parents=True, exist_ok=True)
    
    joined.to_csv(Path(args.out) / "anchors_with_ks.tsv", sep='\t', index=False)
    print(f"Saved: {Path(args.out) / 'anchors_with_ks.tsv'}")

    if 'ks' in joined.columns and joined['ks'].notna().any():
        print("Plotting...")
        plot_anchors_ks(joined, args.out)
        print("Plots saved in:", args.out)
    else:
        print("No Ks values found for anchors; skipping plots.")

if __name__ == '__main__':
    main()
