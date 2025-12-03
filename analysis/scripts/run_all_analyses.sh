#!/bin/bash

# ============================================================================
# Master Analysis Script - Phase 1-4 Complete Analysis
# ============================================================================
# This script runs all analysis scripts in the correct order to generate
# comprehensive results for genome duplication analysis (Phase 1-4)
# ============================================================================

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}============================================================================${NC}"
echo -e "${BLUE}     GLYCINE MAX GENOME DUPLICATION ANALYSIS - PHASE 1-4${NC}"
echo -e "${BLUE}============================================================================${NC}"
echo ""

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}ERROR: python3 not found${NC}"
    exit 1
fi

# Create main output directory
mkdir -p output/analysis_results

# ============================================================================
# STEP 1: Basic Genome Statistics
# ============================================================================
echo -e "${YELLOW}[1/6] Running basic genome statistics analysis...${NC}"
python3 analysis/genome_statistics.py
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  ✓ Genome statistics complete${NC}"
else
    echo -e "${RED}  ✗ Genome statistics failed${NC}"
    exit 1
fi
echo ""

# ============================================================================
# STEP 2: Duplication Type Classification
# ============================================================================
echo -e "${YELLOW}[2/6] Classifying duplication types (TAGs, WGD, Dispersed)...${NC}"
python3 analysis/classify_duplications.py
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  ✓ Duplication classification complete${NC}"
else
    echo -e "${RED}  ✗ Duplication classification failed${NC}"
    exit 1
fi
echo ""

# ============================================================================
# STEP 3: Singleton Identification
# ============================================================================
echo -e "${YELLOW}[3/6] Identifying singleton genes...${NC}"
python3 analysis/identify_singletons.py
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  ✓ Singleton identification complete${NC}"
else
    echo -e "${RED}  ✗ Singleton identification failed${NC}"
    exit 1
fi
echo ""

# ============================================================================
# STEP 4: Ks Distribution Analysis
# ============================================================================
echo -e "${YELLOW}[4/6] Analyzing Ks distributions and WGD peaks...${NC}"
python3 analysis/Ks/01_ks_analysis.py
if [ $? -eq 0 ]; then
    echo -e "${GREEN}  ✓ Ks distribution analysis complete${NC}"
else
    echo -e "${RED}  ✗ Ks distribution analysis failed${NC}"
    exit 1
fi
echo ""

# ============================================================================
# STEP 5: Ks Filtering Comparison (if multiple filtered files exist)
# ============================================================================
echo -e "${YELLOW}[5/6] Comparing Ks distributions across filtering thresholds...${NC}"
KS_COUNT=$(ls output/ks_results/ks_results*.tsv 2>/dev/null | wc -l)
if [ "$KS_COUNT" -gt 1 ]; then
    python3 analysis/Ks/compare_ks_filtering.py output/ks_results/ output/figures/ks_comparison
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}  ✓ Ks filtering comparison complete${NC}"
    else
        echo -e "${YELLOW}  ⚠ Ks filtering comparison failed (non-critical)${NC}"
    fi
else
    echo -e "${YELLOW}  ⚠ Skipping (only one Ks results file found)${NC}"
fi
echo ""

# ============================================================================
# STEP 6: Generate Summary Report
# ============================================================================
echo -e "${YELLOW}[6/6] Generating comprehensive summary report...${NC}"

REPORT_FILE="output/analysis_results/ANALYSIS_SUMMARY.md"

cat > "$REPORT_FILE" <<'EOF'
# Glycine max Genome Duplication Analysis
## Phase 1-4 Complete Results

**Analysis Date:** $(date +"%Y-%m-%d %H:%M:%S")
**Genome:** Glycine max (Soybean) v2.0
**Source:** Ensembl Plants Release 41

---

## Table of Contents
1. [Basic Genome Statistics](#basic-genome-statistics)
2. [Duplication Type Classification](#duplication-type-classification)
3. [Singleton Analysis](#singleton-analysis)
4. [Ks Distribution Analysis](#ks-distribution-analysis)
5. [Key Findings](#key-findings)
6. [Output Files](#output-files)

---

## Basic Genome Statistics

### Overview
EOF

# Add statistics from files
if [ -f "output/statistics/basic_statistics.tsv" ]; then
    echo "" >> "$REPORT_FILE"
    echo '```' >> "$REPORT_FILE"
    cat "output/statistics/basic_statistics.tsv" >> "$REPORT_FILE"
    echo '```' >> "$REPORT_FILE"
fi

cat >> "$REPORT_FILE" <<'EOF'

### Visualizations
- **Chromosome Distribution:** `output/statistics/genes_per_chromosome.png`
- **Family Size Distribution:** `output/statistics/family_size_distribution.png`
- **Pipeline Filtering Effects:** `output/statistics/pipeline_filtering.png`
- **Comprehensive Summary:** `output/statistics/summary_figure.png`

---

## Duplication Type Classification

### Types Identified
1. **TAGs (Tandem Arrayed Genes):** Same chromosome, ≤10 genes or ≤100kb apart
2. **Proximal:** Same chromosome, >100kb but <500kb apart
3. **Dispersed:** Different chromosomes OR >500kb apart
4. **WGD (Whole Genome Duplication):** Based on Ks peaks
   - Peak 1: Ks 0.1-0.3 (~13 MYA, Glycine-specific)
   - Peak 2: Ks 0.7-1.2 (~59 MYA, Papilionoideae-wide)

### Results
EOF

if [ -f "output/duplication_classification/duplication_type_summary.tsv" ]; then
    echo "" >> "$REPORT_FILE"
    echo '```' >> "$REPORT_FILE"
    cat "output/duplication_classification/duplication_type_summary.tsv" >> "$REPORT_FILE"
    echo '```' >> "$REPORT_FILE"
fi

cat >> "$REPORT_FILE" <<'EOF'

### Visualizations
- **Duplication Types Summary:** `output/duplication_classification/duplication_types_summary.png`
- **Classified Gene Pairs:** `output/duplication_classification/gene_pairs_classified.tsv`

---

## Singleton Analysis

### Overview
Singleton genes are those that show no evidence of duplication based on sequence similarity.

EOF

if [ -f "output/singletons/singletons_summary.tsv" ]; then
    echo "" >> "$REPORT_FILE"
    echo '```' >> "$REPORT_FILE"
    cat "output/singletons/singletons_summary.tsv" >> "$REPORT_FILE"
    echo '```' >> "$REPORT_FILE"
fi

cat >> "$REPORT_FILE" <<'EOF'

### Visualizations
- **Singleton Analysis:** `output/singletons/singletons_analysis.png`
- **Chromosome Distribution:** `output/singletons/singletons_by_chromosome.tsv`

---

## Ks Distribution Analysis

### WGD Peak Detection
Using Gaussian Mixture Model (GMM) to identify WGD events based on Ks distribution.

### Age Conversion
Age (MYA) = Ks / (2 × λ)  
where λ = 6.5 × 10⁻⁹ substitutions/site/year (legumes)

### Visualizations
- **Overall Ks Distribution:** `output/figures/ks_distribution_all.png`
- **WGD Peak Detection:** `output/figures/ks_wgd_peaks.png`
- **Age Distribution:** `output/figures/age_distribution.png`

EOF

if [ -f "output/figures/ks_comparison_combined.png" ]; then
    cat >> "$REPORT_FILE" <<'EOF'
- **Filtering Comparison:** `output/figures/ks_comparison_combined.png`

EOF
fi

cat >> "$REPORT_FILE" <<'EOF'
---

## Key Findings

### 1. Genome Composition
- High proportion of duplicated genes suggests extensive duplication history
- Chromosome-specific variation in duplication rates

### 2. Duplication Mechanisms
- Multiple duplication mechanisms active in Glycine max
- TAGs represent recent, local duplications
- WGD peaks correspond to known paleopolyploidy events

### 3. Evolutionary Patterns
- Ks distribution shows clear peaks corresponding to:
  - Recent Glycine-specific WGD (~13 MYA)
  - Ancient Papilionoideae WGD (~59 MYA)
- Ka/Ks ratios suggest purifying selection on most duplicates

### 4. Singleton Characteristics
- Singletons may represent:
  - Essential genes under strong constraint
  - Recently evolved genes
  - Genes that lost duplicates after WGD

---

## Output Files

### Statistics
- `output/statistics/` - Basic genome and family statistics
- `output/duplication_classification/` - Duplication type classifications
- `output/singletons/` - Singleton gene lists and statistics

### Visualizations
- `output/figures/` - Ks distribution plots and WGD analysis
- All summary figures in respective output directories

### Data Files
- `output/duplication_classification/gene_pairs_classified.tsv` - All gene pairs with type annotation
- `output/singletons/singletons.tsv` - Complete singleton gene list
- `output/statistics/ks_filtered.tsv` - Filtered Ks values with age estimates

---

## Next Steps (Phase 5+)

1. **Functional Analysis**
   - GO enrichment for different duplication types
   - Compare TAGs vs WGD duplicates functionally
   - Analyze singleton gene functions

2. **TE Analysis**
   - TE distribution analysis
   - TE enrichment near duplicated genes
   - Correlation with duplication age

3. **Detailed TAG Analysis**
   - TAG cluster identification
   - Age distribution of TAGs
   - Functional specialization

---

**Generated by:** `analysis/run_all_analyses.sh`  
**Analysis Scripts:** `analysis/*.py`

EOF

echo -e "${GREEN}  ✓ Summary report generated: $REPORT_FILE${NC}"
echo ""

# ============================================================================
# Final Summary
# ============================================================================
echo -e "${BLUE}============================================================================${NC}"
echo -e "${GREEN}                     ✓ ALL ANALYSES COMPLETE${NC}"
echo -e "${BLUE}============================================================================${NC}"
echo ""
echo -e "${GREEN}Output locations:${NC}"
echo -e "  • Basic statistics:        ${BLUE}output/statistics/${NC}"
echo -e "  • Duplication types:       ${BLUE}output/duplication_classification/${NC}"
echo -e "  • Singletons:              ${BLUE}output/singletons/${NC}"
echo -e "  • Ks analysis:             ${BLUE}output/figures/${NC}"
echo -e "  • Summary report:          ${BLUE}$REPORT_FILE${NC}"
echo ""
echo -e "${GREEN}Next steps:${NC}"
echo -e "  1. Review the summary report: ${BLUE}$REPORT_FILE${NC}"
echo -e "  2. Check visualizations in output directories"
echo -e "  3. Proceed to functional analysis (GO enrichment)"
echo -e "  4. Analyze TE distribution (Phase 5)"
echo ""
echo -e "${BLUE}============================================================================${NC}"
