#!/usr/bin/env python3
"""
Filter MCScanX .collinearity file for significant syntenic blocks and infer syntenic chromosomes.

Usage:
  python analysis/scripts/filter_collinearity.py --in output/mcscanx/glycine.collinearity --out output/mcscanx/glycine.filtered.collinearity \
      --min-score 50 --min-pairs 3 --max-evalue 1e-5

Notes:
- Keeps only blocks with score >= --min-score, at least --min-pairs anchors, and optional max e-value.
- Preserves original file header (lines before first "## Alignment").
- Infers syntenic chromosomes per block by parsing gene names (supports patterns like Glyma01g, chr01.gene, scaffold01|gene).
- Emits a summary comment line per block: "# Syntenic chromosomes: chrX, chrY, ...".
"""
import argparse
import re
from pathlib import Path

HEADER_RE = re.compile(r"^##\s*Alignment\s+(\d+):.*?score=\s*([0-9\.]+).*?e_value=\s*([eE0-9\-\.]+).*?N=\s*(\d+)")
ANCHOR_RE = re.compile(r"^\s*\d+-\s*\d+:\s+(\S+)\s+(\S+)\s+([eE0-9\-\.]+)")

CHR_PATTERNS = [
    re.compile(r"^(?:Glyma|Gm)(\d{1,2})", re.IGNORECASE),
    re.compile(r"^(?:chr|Chr)(\d{1,2})", re.IGNORECASE),
    re.compile(r"^(\d{1,2})[A-Za-z]"),
    re.compile(r"^(scaffold|Scaffold)(\d{1,3})"),
]

def extract_chr(gene: str):
    for pat in CHR_PATTERNS:
        m = pat.match(gene)
        if m:
            return str(m.groups()[-1])
    m = re.match(r"^(\d{1,2})", gene)
    return m.group(1) if m else None


def filter_collinearity(in_path: Path, out_path: Path, min_score: float, min_pairs: int, max_evalue: float | None):
    keep_block = False
    block_chr_set = set()
    wrote_header = False
    header_buffer = []

    with in_path.open('r') as fin, out_path.open('w') as fout:
        for line in fin:
            # Buffer header lines until first alignment header
            if not wrote_header and not line.startswith("## Alignment"):
                header_buffer.append(line)
                continue

            # On first alignment, write buffered header
            if not wrote_header and line.startswith("## Alignment"):
                if header_buffer:
                    fout.writelines(header_buffer)
                wrote_header = True

            # Handle alignment header
            if line.startswith("## Alignment"):
                m = HEADER_RE.match(line)
                keep_block = False
                block_chr_set.clear()
                block_evalue = None
                if m:
                    score = float(m.group(2))
                    try:
                        block_evalue = float(m.group(3))
                    except Exception:
                        block_evalue = None
                    n_pairs = int(m.group(4))
                    if score >= min_score and n_pairs >= min_pairs:
                        if max_evalue is None or (block_evalue is not None and block_evalue <= max_evalue):
                            keep_block = True
                if keep_block:
                    fout.write(line)
                continue

            # Anchor lines
            am = ANCHOR_RE.match(line)
            if am:
                if keep_block:
                    g1, g2 = am.group(1), am.group(2)
                    c1 = extract_chr(g1)
                    c2 = extract_chr(g2)
                    if c1:
                        block_chr_set.add(c1)
                    if c2:
                        block_chr_set.add(c2)
                    fout.write(line)
                continue

            # Blank line denotes end of block
            if line.strip() == "":
                if keep_block:
                    chrs = sorted(block_chr_set)
                    if chrs:
                        fout.write(f"# Syntenic chromosomes: {', '.join(chrs)}\n")
                keep_block = False
                block_chr_set.clear()
                fout.write(line)
                continue

            # Other lines
            if keep_block:
                fout.write(line)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Filter MCScanX .collinearity and infer syntenic chromosomes")
    ap.add_argument("--in", required=True, dest="infile", help="Input .collinearity file")
    ap.add_argument("--out", required=True, dest="outfile", help="Output filtered .collinearity file")
    ap.add_argument("--min-score", type=float, default=50.0, help="Minimum block score to keep (default: 50)")
    ap.add_argument("--min-pairs", type=int, default=3, help="Minimum anchor pairs N per block (default: 3)")
    ap.add_argument("--max-evalue", type=float, default=None, help="Maximum block e_value (optional)")
    args = ap.parse_args()

    in_path = Path(args.infile)
    out_path = Path(args.outfile)

    filter_collinearity(in_path, out_path, args.min_score, args.min_pairs, args.max_evalue)
    print(f"Filtered collinearity written to {out_path}")
