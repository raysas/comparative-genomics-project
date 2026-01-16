#!/usr/bin/env python3
"""
Compare Tandem Array Gene (TAG) findings from multiple approaches and plot overlap.

Inputs (defaults target the workspace paths):
  --distance     Path to TAGs derived by bp distance threshold (TSV with Gene1, Gene2)
  --genecount    Path to TAGs derived by gene count (TSV with Gene1, Gene2)
  --mcscanx      Path to MCScanX tandem file (CSV pairs per line)
    --tagpairs     Optional path to `TAG_gene_pairs.tsv` (columns geneA, geneB, ...)
    --pairwise     If set, compare only `--tagpairs` vs `--genecount` (2-set overlap)
  --outdir       Output root directory (stats/ and plots/ are created under this)
  --prefix       Output filename prefix (default: tags_overlap)

Outputs:
  - stats/tags_overlap_stats.tsv: set sizes and intersections (pairwise + triple)
  - stats/tags_overlap_membership.tsv: per-gene membership across approaches
    - plots/tags_overlap_venn.png: Venn diagram for 3 sets (if matplotlib-venn available)
  - plots/tags_overlap_counts.png: Fallback counts bar chart if venn not available
    - plots/tags_pairwise_venn.png: Venn diagram for 2 sets (pairwise mode)
    - plots/tags_pairwise_counts.png: Fallback counts bar chart for 2 sets (pairwise mode)

Approach:
  - Build gene sets from each source by flattening all gene pairs.
  - Compute set sizes, pairwise intersections, and triple intersection.
  - Produce overlap visualization.
"""
import argparse
import csv
from pathlib import Path
from typing import Set, Tuple, Dict

import pandas as pd
import matplotlib.pyplot as plt


def read_pairs_tsv(path: Path, gene1_col: str = "Gene1", gene2_col: str = "Gene2") -> Set[str]:
    """Read a TSV of gene pairs and return the set of unique genes.

    Expected columns include gene identifiers in `gene1_col` and `gene2_col`.
    Falls back to case-insensitive matching if needed.
    """
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    g1 = gene1_col.lower()
    g2 = gene2_col.lower()
    if g1 not in cols or g2 not in cols:
        # try alternatives
        candidates1 = ["gene1", "peptide1", "id1", "query"]
        candidates2 = ["gene2", "peptide2", "id2", "subject"]
        for cand in candidates1:
            if cand in cols:
                gene1_col = cols[cand]
                break
        for cand in candidates2:
            if cand in cols:
                gene2_col = cols[cand]
                break
    genes = set(df[gene1_col].astype(str)) | set(df[gene2_col].astype(str))
    return genes


def read_mcscanx_tandem(path: Path) -> Set[str]:
    """Read MCScanX tandem file where each line has two gene IDs separated by comma.
    Returns the set of unique genes participating in tandem pairs.
    """
    genes: Set[str] = set()
    with path.open("r", newline="") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if not row:
                continue
            # Handle potential whitespace
            items = [item.strip() for item in (row if len(row) > 1 else row[0].split(","))]
            if len(items) >= 2:
                genes.add(items[0])
                genes.add(items[1])
    return genes


def read_tag_gene_pairs(path: Path) -> Set[str]:
    """Read TAG_gene_pairs.tsv with columns geneA, geneB, tag_id, etc., return unique genes."""
    df = pd.read_csv(path, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    g1 = cols.get("genea", None)
    g2 = cols.get("geneb", None)
    if not g1 or not g2:
        raise ValueError("TAG_gene_pairs.tsv must contain columns 'geneA' and 'geneB'")
    return set(df[g1].astype(str)) | set(df[g2].astype(str))


def compute_overlap_sets(sets: Dict[str, Set[str]]) -> Dict[str, int]:
    """Compute sizes and intersections for three sets.

    Returns a dict with counts for each set and intersections.
    Keys:
      - size_<label>
      - inter_<labelA>_<labelB>
      - inter_all
    """
    labels = list(sets.keys())
    if len(labels) != 3:
        raise ValueError("Expected exactly three sets for overlap computation")
    A, B, C = sets[labels[0]], sets[labels[1]], sets[labels[2]]
    counts = {
        f"size_{labels[0]}": len(A),
        f"size_{labels[1]}": len(B),
        f"size_{labels[2]}": len(C),
        f"inter_{labels[0]}_{labels[1]}": len(A & B),
        f"inter_{labels[0]}_{labels[2]}": len(A & C),
        f"inter_{labels[1]}_{labels[2]}": len(B & C),
        "inter_all": len(A & B & C),
    }
    return counts


def save_membership_table(sets: Dict[str, Set[str]], out_path: Path) -> None:
    """Save per-gene membership across approaches to TSV."""
    labels = list(sets.keys())
    all_genes = sorted(set().union(*sets.values()))
    rows = []
    for g in all_genes:
        row = {"Gene": g}
        present = 0
        for lbl in labels:
            in_set = g in sets[lbl]
            row[lbl] = in_set
            present += int(in_set)
        row["Overlap_Count"] = present
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def plot_venn_or_counts(sets: Dict[str, Set[str]], plots_dir: Path, prefix: str) -> Tuple[Path, Path]:
    """Plot Venn diagram (if matplotlib-venn available) else fallback counts bar chart.

    Returns tuple of (venn_path or None, counts_path or None).
    """
    labels = list(sets.keys())
    A, B, C = sets[labels[0]], sets[labels[1]], sets[labels[2]]
    venn_path = plots_dir / f"{prefix}_venn.png"
    counts_path = plots_dir / f"{prefix}_counts.png"

    try:
        from matplotlib_venn import venn3
        fig, ax = plt.subplots(figsize=(8, 6))
        venn3(subsets=(len(A - (B | C)),
                       len(B - (A | C)),
                       len(A & B - C),
                       len(C - (A | B)),
                       len(A & C - B),
                       len(B & C - A),
                       len(A & B & C)),
              set_labels=labels)
        ax.set_title("Overlap of TAG approaches")
        fig.tight_layout()
        fig.savefig(venn_path, dpi=200)
        plt.close(fig)
        return venn_path, None
    except Exception:
        # Fallback: simple counts bar chart
        sizes = [len(A), len(B), len(C), len(A & B), len(A & C), len(B & C), len(A & B & C)]
        names = [f"{labels[0]}", f"{labels[1]}", f"{labels[2]}",
                 f"{labels[0]}∩{labels[1]}", f"{labels[0]}∩{labels[2]}", f"{labels[1]}∩{labels[2]}",
                 "all three"]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(range(len(sizes)), sizes, color=["#6baed6", "#9ecae1", "#c6dbef", "#74c476", "#31a354", "#006d2c", "#756bb1"])
        ax.set_xticks(range(len(sizes)))
        ax.set_xticklabels(names, rotation=30, ha="right")
        ax.set_ylabel("Gene count")
        ax.set_title("TAG approach counts and overlaps")
        fig.tight_layout()
        fig.savefig(counts_path, dpi=200)
        plt.close(fig)
        return None, counts_path


def plot_pairwise(A: Set[str], B: Set[str], labels: Tuple[str, str], plots_dir: Path, prefix: str) -> Tuple[Path, Path]:
    """Plot 2-set Venn if available, else fallback bar chart."""
    venn_path = plots_dir / f"{prefix}_pairwise_venn.png"
    counts_path = plots_dir / f"{prefix}_pairwise_counts.png"
    try:
        from matplotlib_venn import venn2
        fig, ax = plt.subplots(figsize=(7, 5))
        venn2(subsets=(len(A - B), len(B - A), len(A & B)), set_labels=list(labels))
        ax.set_title("Pairwise overlap")
        fig.tight_layout()
        fig.savefig(venn_path, dpi=200)
        plt.close(fig)
        return venn_path, None
    except Exception:
        sizes = [len(A), len(B), len(A & B)]
        names = [labels[0], labels[1], f"{labels[0]}∩{labels[1]}"]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(range(len(sizes)), sizes, color=["#6baed6", "#9ecae1", "#74c476"])
        ax.set_xticks(range(len(sizes)))
        ax.set_xticklabels(names, rotation=0, ha="center")
        ax.set_ylabel("Gene count")
        ax.set_title("Pairwise TAG counts and overlap")
        fig.tight_layout()
        fig.savefig(counts_path, dpi=200)
        plt.close(fig)
        return None, counts_path


def main():
    parser = argparse.ArgumentParser(description="Compare TAG approaches and plot overlaps")
    parser.add_argument("--distance", type=Path, default=Path("analysis/duplication_types/TAGs/TAGs_distance_100000bp.tsv"),
                        help="TSV with Gene1/Gene2 for distance-based TAGs")
    parser.add_argument("--genecount", type=Path, default=Path("analysis/duplication_types/TAGs/TAGs_genecount_10genes.tsv"),
                        help="TSV with Gene1/Gene2 for gene-count-based TAGs")
    parser.add_argument("--mcscanx", type=Path, default=Path("analysis/duplication_types/TAGs/glycine.tandem"),
                        help="MCScanX tandem pairs file (CSV per line)")
    parser.add_argument("--tagpairs", type=Path, default=None,
                        help="TAG_gene_pairs.tsv path (with geneA/geneB)")
    parser.add_argument("--pairwise", action="store_true",
                        help="Compare only TAG_gene_pairs vs genecount (2 sets)")
    parser.add_argument("--outdir", type=Path, default=Path("analysis/duplication_types/TAGs"),
                        help="Output root dir; stats/ and plots/ will be created")
    parser.add_argument("--prefix", default="tags_overlap", help="Output file prefix")
    args = parser.parse_args()

    # Prepare outputs
    stats_dir = args.outdir / "stats"
    plots_dir = args.outdir / "plots"
    stats_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    if args.pairwise:
        if args.tagpairs is None:
            raise SystemExit("--pairwise requires --tagpairs to be provided")
        genecount_genes = read_pairs_tsv(args.genecount)
        tagpairs_genes = read_tag_gene_pairs(args.tagpairs)
        # Stats
        rows = [{
            "size_tagpairs": len(tagpairs_genes),
            "size_genecount": len(genecount_genes),
            "inter_tagpairs_genecount": len(tagpairs_genes & genecount_genes),
        }]
        pd.DataFrame(rows).to_csv(stats_dir / f"{args.prefix}_pairwise_stats.tsv", sep="\t", index=False)
        # Membership table
        sets_pw = {"tagpairs": tagpairs_genes, "genecount": genecount_genes}
        save_membership_table(sets_pw, stats_dir / f"{args.prefix}_pairwise_membership.tsv")
        # Plot
        venn_path, counts_path = plot_pairwise(tagpairs_genes, genecount_genes, ("TAG_gene_pairs", "genecount"), plots_dir, args.prefix)
        print("Saved stats:", stats_dir / f"{args.prefix}_pairwise_stats.tsv")
        print("Saved membership:", stats_dir / f"{args.prefix}_pairwise_membership.tsv")
        if venn_path:
            print("Saved Venn:", venn_path)
        if counts_path:
            print("Saved counts chart:", counts_path)
    else:
        # Read sets
        distance_genes = read_pairs_tsv(args.distance)
        genecount_genes = read_pairs_tsv(args.genecount)
        mcscanx_genes = read_mcscanx_tandem(args.mcscanx)

        sets = {
            "distance": distance_genes,
            "genecount": genecount_genes,
            "mcscanx": mcscanx_genes,
        }

        # Compute and save stats
        counts = compute_overlap_sets(sets)
        pd.DataFrame([counts]).to_csv(stats_dir / f"{args.prefix}_stats.tsv", sep="\t", index=False)

        # Save membership table
        save_membership_table(sets, stats_dir / f"{args.prefix}_membership.tsv")

        # Plot
        venn_path, counts_path = plot_venn_or_counts(sets, plots_dir, args.prefix)
        print("Saved stats:", stats_dir / f"{args.prefix}_stats.tsv")
        print("Saved membership:", stats_dir / f"{args.prefix}_membership.tsv")
        if venn_path:
            print("Saved Venn:", venn_path)
        if counts_path:
            print("Saved counts chart:", counts_path)


if __name__ == "__main__":
    main()
