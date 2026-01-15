#!/usr/bin/env python3
"""
Analyze syntenic chromosome relationships and assign anchors to the two WGD events using Ks.

Inputs:
  --anchors  TSV from MCScanX anchors joined with Ks (e.g., analysis/mcscanx_anchors/anchors_with_ks.tsv)
             Expected columns: block,gene1,gene2,ks (others optional)
  --metadata Protein metadata with genomic positions (default: data/glycine_max/processed/protein_info.csv)
             Must include: peptide_id,chromosome,start_pos,end_pos
  --outdir   Output directory for results

WGD assignment options:
  - Default: Fit a 2-component GaussianMixture on Ks and assign each pair to the highest responsibility component.
  - Optional thresholds: use --wgd1-max and --wgd2-min to classify by Ks ranges instead of mixture.

Outputs:
  - chromosome_pairs_wgd.tsv: chr1,chr2,WGD,label counts, mean/median Ks, responsibility proportion
  - chromosome_summary_wgd.tsv: per chr1, top matching chr2 per WGD with counts and median Ks
  - wgd_mixture_stats.tsv: means, sds, weights (if mixture used)
  - Optional heatmaps (PNG): anchors_by_chr_wgd.png
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
import os


def load_anchors(anchors_path: Path) -> pd.DataFrame:
    df = pd.read_csv(anchors_path, sep='\t')
    cols = {c.lower(): c for c in df.columns}
    required = ['gene1','gene2']
    for r in required:
        if r not in cols:
            raise ValueError(f"Anchors file missing required column: {r}")
    # Normalize column names of interest
    out = df.rename(columns={
        cols['gene1']: 'gene1',
        cols['gene2']: 'gene2',
    })
    if 'block' in cols:
        out = out.rename(columns={cols['block']: 'block'})
    if 'ks' in cols:
        out = out.rename(columns={cols['ks']: 'ks'})
    else:
        raise ValueError("Anchors file missing 'ks' column for WGD assignment")
    return out


def load_metadata(meta_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(meta_path)
    colmap = {c.lower(): c for c in meta.columns}
    required = ['peptide_id','chromosome','start_pos','end_pos']
    for r in required:
        if r not in colmap:
            raise ValueError(f"Metadata missing required column: {r}")
    meta = meta.rename(columns={
        colmap['peptide_id']: 'peptide_id',
        colmap['chromosome']: 'chromosome',
        colmap['start_pos']: 'start_pos',
        colmap['end_pos']: 'end_pos',
    })
    meta['peptide_id'] = meta['peptide_id'].astype(str)
    meta['mid_pos'] = (meta['start_pos'].astype(float) + meta['end_pos'].astype(float)) / 2.0
    return meta[['peptide_id','chromosome','mid_pos']]


def map_genes_to_chromosomes(anchors: pd.DataFrame, meta: pd.DataFrame, chr_format: str | None) -> pd.DataFrame:
    idx = meta.set_index('peptide_id')
    rows = []
    for _, r in anchors.iterrows():
        g1 = str(r['gene1']) if pd.notna(r['gene1']) else None
        g2 = str(r['gene2']) if pd.notna(r['gene2']) else None
        if g1 is None or g2 is None:
            continue
        if g1 not in idx.index or g2 not in idx.index:
            continue
        chr1 = idx.at[g1, 'chromosome']
        chr2 = idx.at[g2, 'chromosome']
        # Optional normalized chromosome formatting
        def fmt(c):
            if chr_format:
                try:
                    return chr_format % int(c)
                except Exception:
                    return str(c)
            return str(c)
        rows.append({
            'chr1': fmt(chr1),
            'chr2': fmt(chr2),
            'ks': float(r['ks']) if pd.notna(r['ks']) else np.nan,
            'block': r['block'] if 'block' in anchors.columns else None,
        })
    out = pd.DataFrame(rows)
    # Filter reasonable Ks
    out = out[out['ks'].between(0, 2, inclusive='both')]
    return out


def assign_wgd_by_mixture(ks_values: np.ndarray) -> dict:
    ks = ks_values.reshape(-1, 1)
    gmm = GaussianMixture(n_components=2, random_state=42, max_iter=1000)
    gmm.fit(ks)
    means = gmm.means_.flatten()
    variances = gmm.covariances_.flatten()
    weights = gmm.weights_.flatten()
    order = np.argsort(means)
    means = means[order]
    sds = np.sqrt(variances[order])
    weights = weights[order]
    resp = gmm.predict_proba(ks)[:, order]  # columns align to sorted components
    labels_idx = resp.argmax(axis=1)
    labels = np.array(['WGD1_lowKs','WGD2_highKs'])[labels_idx]
    return {
        'labels': labels,
        'responsibilities': resp.max(axis=1),
        'means': means,
        'sds': sds,
        'weights': weights,
    }


def assign_wgd_by_thresholds(ks_values: np.ndarray, wgd1_max: float, wgd2_min: float) -> dict:
    ks = ks_values
    labels = np.where(ks <= wgd1_max, 'WGD1_lowKs', np.where(ks >= wgd2_min, 'WGD2_highKs', 'Unassigned'))
    resp = np.ones_like(ks, dtype=float)
    return {
        'labels': labels,
        'responsibilities': resp,
        'means': np.array([np.nan, np.nan]),
        'sds': np.array([np.nan, np.nan]),
        'weights': np.array([np.nan, np.nan]),
    }


def summarize_pairs(df: pd.DataFrame) -> pd.DataFrame:
    grp = df.groupby(['chr1','chr2','wgd_label'])
    out = grp.agg(
        n_anchors=('ks','count'),
        mean_ks=('ks','mean'),
        median_ks=('ks','median'),
        mean_resp=('responsibility','mean')
    ).reset_index()
    # Also add symmetric pairs (optional): ensure consistent ordering
    return out


def summarize_per_chr(df: pd.DataFrame) -> pd.DataFrame:
    grp = df.groupby(['chr1','wgd_label','chr2']).agg(n=('ks','count'), med=('ks','median')).reset_index()
    # For each (chr1, wgd), pick top chr2 by count
    tops = grp.sort_values(['chr1','wgd_label','n','med'], ascending=[True, True, False, True]).groupby(['chr1','wgd_label']).head(1)
    return tops.rename(columns={'n':'top_count','med':'top_median_ks'})


def main():
    ap = argparse.ArgumentParser(description='Analyze syntenic chromosome matches and assign to two WGDs by Ks')
    ap.add_argument('--anchors', required=True, help='Anchors with Ks TSV')
    ap.add_argument('--metadata', default='data/glycine_max/processed/protein_info.csv', help='Protein metadata CSV')
    ap.add_argument('--outdir', required=True, help='Output directory')
    ap.add_argument('--use-thresholds', action='store_true', help='Use Ks thresholds instead of mixture')
    ap.add_argument('--wgd1-max', type=float, default=0.3, help='Max Ks for WGD1 (low Ks)')
    ap.add_argument('--wgd2-min', type=float, default=0.45, help='Min Ks for WGD2 (high Ks)')
    ap.add_argument('--chr-format', type=str, default='chr%02d', help='Printf-style format for chromosome names (e.g., chr%02d)')
    args = ap.parse_args()

    anchors_path = Path(args.anchors)
    meta_path = Path(args.metadata)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    anchors = load_anchors(anchors_path)
    meta = load_metadata(meta_path)
    mapped = map_genes_to_chromosomes(anchors, meta, chr_format=args.chr_format if args.chr_format else None)
    if mapped.empty:
        raise RuntimeError('No mapped anchors with Ks found after filtering')

    ks = mapped['ks'].to_numpy()
    if args.use_thresholds:
        assign = assign_wgd_by_thresholds(ks, args.wgd1_max, args.wgd2_min)
        mixture_stats = pd.DataFrame({'metric':['means','sds','weights'], 'comp1':[np.nan,np.nan,np.nan], 'comp2':[np.nan,np.nan,np.nan]})
    else:
        assign = assign_wgd_by_mixture(ks)
        mixture_stats = pd.DataFrame({
            'component': ['WGD1_lowKs','WGD2_highKs'],
            'mean': assign['means'],
            'sd': assign['sds'],
            'weight': assign['weights'],
        })
    mapped['wgd_label'] = assign['labels']
    mapped['responsibility'] = assign['responsibilities']

    pairs_summary = summarize_pairs(mapped)
    per_chr_summary = summarize_per_chr(mapped)

    pairs_summary.to_csv(outdir / 'chromosome_pairs_wgd.tsv', sep='\t', index=False)
    per_chr_summary.to_csv(outdir / 'chromosome_summary_wgd.tsv', sep='\t', index=False)
    mixture_stats.to_csv(outdir / 'wgd_mixture_stats.tsv', sep='\t', index=False)

    print(f"Wrote: {outdir / 'chromosome_pairs_wgd.tsv'}")
    print(f"Wrote: {outdir / 'chromosome_summary_wgd.tsv'}")
    print(f"Wrote: {outdir / 'wgd_mixture_stats.tsv'}")


if __name__ == '__main__':
    main()
