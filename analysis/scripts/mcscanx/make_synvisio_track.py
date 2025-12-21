#!/usr/bin/env python3
"""
Make SynVisio track files from anchors_with_ks.tsv

Supports two track types:
- binned median-Ks per chromosome window (recommended)
- binned counts (anchor density)

Track formats written (TSV):
- Optional first line: "min=0,max=1" to force SynVisio scaling
- Then rows: "chr\tstart\tend\tvalue"

Usage examples:
  python analysis/scripts/make_synvisio_track.py \
    --in analysis/mcscanx_anchors/anchors_with_ks.tsv \
    --out output/mcscanx/tracks/anchors_medianKs_100kb.tsv \
    --mode median-ks --bin-size 100000 --min 0 --max 1

  python analysis/scripts/make_synvisio_track.py \
    --in analysis/mcscanx_anchors/anchors_with_ks.tsv \
    --out output/mcscanx/tracks/anchors_density_100kb.tsv \
    --mode count --bin-size 100000
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

CHR_CANDIDATES = [
    'chr','chrom','chromosome',
    'gene_chr','gene1_chr','gene2_chr',
    'anchor_chr'
]
POS_CANDIDATES = [
    'pos','position','start','coord',
    'gene_pos','gene1_pos','gene2_pos',
    'anchor_pos','midpoint'
]


def find_cols(df: pd.DataFrame, candidates: list[str]) -> list[str]:
    cols = []
    lower_map = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c in lower_map:
            cols.append(lower_map[c])
    # also try exact names as-is
    for c in candidates:
        if c in df.columns and c not in cols:
            cols.append(c)
    return cols


def load_metadata(meta_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(meta_path)
    # Expected columns: peptide_id,chromosome,start_pos,end_pos
    # Allow flexible capitalization
    cols = {c.lower(): c for c in meta.columns}
    required = ['peptide_id','chromosome','start_pos','end_pos']
    for r in required:
        if r not in cols:
            raise ValueError(f"Metadata missing required column: {r}")
    meta = meta.rename(columns={
        cols['peptide_id']:'peptide_id',
        cols['chromosome']:'chromosome',
        cols['start_pos']:'start_pos',
        cols['end_pos']:'end_pos',
    })
    # Compute midpoint position
    meta['mid_pos'] = ((meta['start_pos'].astype(float) + meta['end_pos'].astype(float)) / 2.0)
    # Ensure peptide_id is string for joins
    meta['peptide_id'] = meta['peptide_id'].astype(str)
    return meta[['peptide_id','chromosome','mid_pos']]


def longify_coords(df: pd.DataFrame, metadata: pd.DataFrame | None = None, chr_format: str | None = None) -> pd.DataFrame:
    # Try to build a long format of (chr, pos, ks)
    cols_chr = find_cols(df, CHR_CANDIDATES)
    cols_pos = find_cols(df, POS_CANDIDATES)
    # Case 1: direct columns present
    if cols_chr and cols_pos:
        pass
    else:
        # Case 2: derive from gene columns using metadata
        if metadata is None:
            raise ValueError("Could not detect chromosome/position columns and no metadata provided. Pass --metadata protein_info.csv")
        # gene columns to use
        gene_cols = [c for c in df.columns if c.lower() in ('gene','gene1','gene2')]
        if not gene_cols:
            raise ValueError("Input file lacks gene columns (gene1/gene2) needed to map to coordinates.")
        # Build long from genes
        records = []
        # Find ks column if any
        ks_col = None
        for c in df.columns:
            if c.lower() == 'ks':
                ks_col = c
                break
        meta_idx = metadata.set_index('peptide_id')
        for _, r in df.iterrows():
            ks_val = float(r[ks_col]) if (ks_col and pd.notna(r[ks_col])) else np.nan
            for gc in gene_cols:
                gid = r[gc]
                if pd.isna(gid):
                    continue
                gid = str(gid)
                if gid in meta_idx.index:
                    chrom = meta_idx.at[gid, 'chromosome']
                    pos = float(meta_idx.at[gid, 'mid_pos'])
                    # Optionally format chromosome name
                    if chr_format:
                        try:
                            chrom_fmt = chr_format % int(chrom)
                        except Exception:
                            chrom_fmt = str(chrom)
                    else:
                        chrom_fmt = str(chrom)
                    rec = {'chr': chrom_fmt, 'pos': pos}
                    if not np.isnan(ks_val):
                        rec['ks'] = ks_val
                    records.append(rec)
        if not records:
            raise ValueError("No coordinates could be derived from metadata for the provided genes.")
        return pd.DataFrame.from_records(records)

    # collect pairs of chr/pos columns that correspond to the same gene index
    pairs = []
    for chrc in cols_chr:
        base = chrc.replace('_chr','').replace('chromosome','').replace('chrom','').replace('chr','')
        # try find matching pos
        candidates = [p for p in cols_pos if base and base in p]
        if not candidates:
            # fallback: if only one pos column exists, use it
            if len(cols_pos) == 1:
                candidates = [cols_pos[0]]
        for pc in candidates:
            pairs.append((chrc, pc))
    if not pairs:
        # last resort: use first chr and first pos
        pairs = [(cols_chr[0], cols_pos[0])]

    records = []
    has_ks = 'ks' in {c.lower(): c for c in df.columns}
    ks_col = None
    if has_ks:
        for c in df.columns:
            if c.lower() == 'ks':
                ks_col = c
                break
    for _, r in df.iterrows():
        for chrc, pc in pairs:
            chrval = r.get(chrc)
            posval = r.get(pc)
            if pd.notna(chrval) and pd.notna(posval):
                try:
                    posf = float(posval)
                except Exception:
                    continue
                rec = {
                    'chr': str(chrval),
                    'pos': posf,
                }
                if ks_col is not None and pd.notna(r[ks_col]):
                    rec['ks'] = float(r[ks_col])
                records.append(rec)
    if not records:
        raise ValueError("No usable (chr,pos) records found after parsing.")
    out = pd.DataFrame.from_records(records)
    # Normalize chr names to a common string (e.g., leading 'chr') only if consistent – otherwise leave as-is
    return out


def write_track_interval(df_bins: pd.DataFrame, out_path: Path, min_val: float|None, max_val: float|None):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w') as w:
        if min_val is not None and max_val is not None:
            w.write(f"min={min_val},max={max_val}\n")
        for _, r in df_bins.iterrows():
            w.write(f"{r['chr']}\t{int(r['start'])}\t{int(r['end'])}\t{float(r['value'])}\n")


def bin_and_summarize(df_long: pd.DataFrame, bin_size: int, mode: str) -> pd.DataFrame:
    rows = []
    for chrom, g in df_long.groupby('chr'):
        if g.empty:
            continue
        max_pos = int(np.nanmax(g['pos']))
        if max_pos < 1:
            continue
        bins = np.arange(0, max_pos + bin_size, bin_size)
        if len(bins) < 2:
            bins = np.array([0, max(1, bin_size)])
        idx = np.digitize(g['pos'].to_numpy(), bins, right=True)
        for i in range(1, len(bins)):
            start = int(bins[i-1]) + 1
            end = int(bins[i])
            sub = g[idx == i]
            if mode == 'count':
                val = float(len(sub))
            elif mode == 'median-ks':
                if 'ks' in sub.columns and not sub['ks'].empty:
                    val = float(sub['ks'].median()) if sub['ks'].notna().any() else np.nan
                else:
                    val = np.nan
            else:
                raise ValueError("Unknown mode; use 'count' or 'median-ks'")
            if not np.isnan(val):
                rows.append({'chr': chrom, 'start': start, 'end': end, 'value': val})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description='Create SynVisio track from anchors_with_ks.tsv')
    ap.add_argument('--in', required=True, dest='infile', help='Input TSV (anchors_with_ks.tsv)')
    ap.add_argument('--out', required=True, dest='outfile', help='Output track TSV path')
    ap.add_argument('--mode', choices=['median-ks','count'], default='median-ks', help='Summary per bin')
    ap.add_argument('--bin-size', type=int, default=100_000, help='Bin size in bp (default 100kb)')
    ap.add_argument('--min', type=float, default=None, help='Optional SynVisio scale min (first-line header)')
    ap.add_argument('--max', type=float, default=None, help='Optional SynVisio scale max (first-line header)')
    ap.add_argument('--metadata', type=str, default='data/glycine_max/processed/protein_info.csv', help='Protein metadata with coordinates (peptide_id, chromosome, start_pos, end_pos)')
    ap.add_argument('--chr-format', type=str, default='chr%02d', help='Printf-style chromosome naming, e.g., chr%02d (set to empty to keep as-is)')
    args = ap.parse_args()

    in_path = Path(args.infile)
    out_path = Path(args.outfile)

    df = pd.read_csv(in_path, sep='\t')
    meta_df = None
    if args.metadata:
        meta_path = Path(args.metadata)
        if not meta_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {meta_path}")
        meta_df = load_metadata(meta_path)
    chr_fmt = args.chr_format if args.chr_format else None
    df_long = longify_coords(df, metadata=meta_df, chr_format=chr_fmt)

    # If mode is median-ks but 'ks' is missing, fall back to count
    mode = args.mode
    if mode == 'median-ks' and 'ks' not in df_long.columns:
        print("Warning: 'ks' column not found; falling back to count mode")
        mode = 'count'

    df_bins = bin_and_summarize(df_long, args.bin_size, mode)

    write_track_interval(df_bins, out_path, args.min, args.max)
    print(f"Wrote track: {out_path} ({len(df_bins)} intervals) mode={mode} bin={args.bin_size}")


if __name__ == '__main__':
    main()
