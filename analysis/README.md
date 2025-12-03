
# Analysis Workflows - Glycine Max Comparative Genomics

This directory contains all analysis scripts for the Glycine max comparative genomics project, including genome statistics, duplication classification, Ks-based WGD detection, and MCScanX synteny analysis.


## Directory Structure

```
analysis/
├── Ks/
│   ├── figures/                  # Ks plots and figures
│   ├── ks_comparison/            # Ks filtering comparison results
│   ├── ks_plots/                 # Distribution and scatter plots
│   ├── statistics/               # Filtered Ks and WGD peak tables
├── duplicated_genes/
│   ├── TAGs/
│       ├── 01_identify_TAGs.py   # Identify tandemly arrayed genes (TAGs)
│       ├── 03_TAG_age.R          # TAG age analysis (R)
│   
├── duplication_classification/
│   ├── duplication_types_summary.pdf/png/tsv
│   ├── gene_pairs_classified.tsv
├── mcscanx_anchors/
│   ├── anchors_with_ks.tsv       # MCScanX anchors annotated with Ks
│   ├── block_median_ks.pdf/png   # Median Ks per block
│   ├── ks_hist.pdf/png           # Ks histogram for anchors
├── singletons/
│   ├── singletons.tsv            # Singleton gene list
│   ├── singletons_analysis.pdf/png
│   ├── singletons_by_chromosome.tsv
│   ├── singletons_summary.tsv
├── statistics/
│   ├── basic_statistics.tsv
│   ├── family_size_distribution.pdf/png
│   ├── genes_per_chromosome.pdf/png/tsv
│   ├── pipeline_filtering.pdf/png
│   ├── summary_figure.pdf/png
│   ├── top_20_families.tsv
├── blast_analysis/
│   ├── blast_e-5/
│   ├── blast_summary.txt
│   ├── paralog_counts.tsv
│   ├── top_gene_families.txt
├── scripts/
│   ├── analyze_blast_results.sh
│   ├── annotate_mcscanx_anchors_with_ks.py
│   ├── classify_duplications.py
│   ├── compare_ks_filtering.py
│   ├── genome_statistics.py
│   ├── identify_singletons.py
│   ├── ks_analysis.py
│   ├── plot_ks_results.py
│   ├── plot_ks_wgd.py
│   ├── run_all_analyses.sh
│   ├── run_mcscanx.sh
```

---

## MCScanX Synteny Analysis

**Step 1: Prepare MCScanX inputs**

Use scripts in `scripts/` to convert protein info and BLAST results:

```bash
# Convert protein info to MCScanX GFF format
bash scripts/convert_protein_info_to_gff.sh data/protein_info_longest.csv output/mcscanx/soybean.gff gm

# Filter BLAST for high-quality matches and convert to MCScanX format
bash scripts/convert_blast_for_mcscanx.sh output/blast_output/blast_results_with_coverage.tsv output/mcscanx/soybean.blast
```

**Step 2: Run MCScanX**

```bash
cd output/mcscanx
MCScanX soybean
```

**Step 3: Annotate MCScanX anchors with Ks**

```bash
python3 scripts/annotate_mcscanx_anchors_with_ks.py \
  --col output/mcscanx/soybean.collinearity \
  --ks output/ks_results/ks_results_filtered.tsv \
  --out output/mcscanx/soybean_anchors_with_ks.tsv \
  --plots output/mcscanx/soybean_anchors_plots
```

**Outputs:**
- `soybean.collinearity` - MCScanX syntenic blocks
- `soybean_anchors_with_ks.tsv` - Anchor pairs annotated with Ks
- `ks_hist.png/pdf` - Ks distribution for anchors
- `block_median_ks.png/pdf` - Median Ks per block

---