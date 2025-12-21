#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    Analyzes BLAST results to provide comprehensive statistics
#    before proceeding to coverage computation
# -- Usage:
#    bash ./scripts/analyze_blast_results.sh [-i INPUT_FILE] [-o OUTPUT_DIR]
# -- default (without params):
#    Auto-detects from output/{species}/blast_output/blast_results.tsv
# --------------------------------------------------------------------

# -- default parameters
INPUT_FILE=''
OUTPUT_DIR=''
AUTO_DETECT=true

# -- get arguments
while getopts "i:o:h" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}"; AUTO_DETECT=false ;;
        o) OUTPUT_DIR="${OPTARG}"; AUTO_DETECT=false ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE       Input BLAST results TSV file
  -o DIR        Output directory for analysis reports
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects:
    Input:  output/{species}/blast_output/blast_results.tsv
    Output: output/{species}/blast_output/

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i output/glycine_max/blast_output/blast_results.tsv \\
     -o output/glycine_max/blast_output/

EOF
            exit 0
            ;;
        *)
            echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
    esac
done

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting BLAST results..."
    
    SPECIES_DIRS=(output/*/blast_output)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        OUTPUT_DIR="${SPECIES_DIRS[0]}"
        INPUT_FILE="${OUTPUT_DIR}/blast_results.tsv"
        SPECIES_NAME=$(basename $(dirname "$OUTPUT_DIR"))
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Input:  $INPUT_FILE"
        echo "   Output: $OUTPUT_DIR"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species BLAST results found"
        echo "       Specify file explicitly with -i"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No BLAST results found in output/*/"
        echo "       Run 2_blast.sh first"
        exit 1
    fi
fi

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: BLAST results file not found: $INPUT_FILE"
    exit 1
fi

echo ""
echo "========================================================"
echo "           BLAST Results Analysis"
echo "========================================================"
echo ""

# Output files
## Ensure output directory exists (support cases where -i provided but -o not)
# If OUTPUT_DIR not set, default to the directory containing the input file
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=$(dirname "$INPUT_FILE")
    # if dirname returns '.', convert to current working directory
    [ "$OUTPUT_DIR" = "." ] && OUTPUT_DIR="$PWD"
fi

# Create directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

SUMMARY_FILE="${OUTPUT_DIR}/blast_summary.txt"
PARALOG_COUNTS="${OUTPUT_DIR}/paralog_counts.tsv"
TOP_FAMILIES="${OUTPUT_DIR}/top_gene_families.txt"

# Start summary file
exec > >(tee "$SUMMARY_FILE")

echo "Input file: $INPUT_FILE"
echo "Analysis date: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# ============================================
# Main Analysis AWK Script
# ============================================

awk -v paralog_file="$PARALOG_COUNTS" -v top_file="$TOP_FAMILIES" '
BEGIN {
    print "=== BASIC STATISTICS ==="
    OFS="\t"
}

{
    total_hits++
    queries[$1]++
    targets[$2]++
    
    # Track self-hits vs non-self
    if ($1 == $2) {
        self_hits++
    } else {
        non_self_hits++
        paralog_count[$1]++
    }
    
    # Identity statistics
    identity_sum += $3
    if ($3 >= 90) high_id++
    else if ($3 >= 70) med_id++
    else low_id++
    
    # E-value statistics  
    if ($11 < 1e-50) very_sig++
    else if ($11 < 1e-10) sig++
    else if ($11 < 1e-5) weak++
    else insignificant++
    
    # Alignment length
    aln_len_sum += $4
    if ($4 < 50) short_aln++
    else if ($4 < 200) med_aln++
    else long_aln++
}

END {
    # Print basic stats
    print ""
    printf "Total BLAST hits:          %d\n", total_hits
    printf "  Self-hits:               %d (%.1f%%)\n", self_hits, 100*self_hits/total_hits
    printf "  Non-self hits:           %d (%.1f%%)\n", non_self_hits, 100*non_self_hits/total_hits
    
    print ""
    printf "Unique query proteins:     %d\n", length(queries)
    printf "Unique target proteins:    %d\n", length(targets)
    
    # Identity distribution
    print ""
    print "=== PERCENT IDENTITY DISTRIBUTION ==="
    printf "Average %% identity:        %.2f%%\n", identity_sum/total_hits
    print ""
    printf "  High (≥90%%):             %d (%.1f%%)\n", high_id, 100*high_id/total_hits
    printf "  Medium (70-90%%):         %d (%.1f%%)\n", med_id, 100*med_id/total_hits
    printf "  Low (<70%%):              %d (%.1f%%)\n", low_id, 100*low_id/total_hits
    
    # E-value distribution
    print ""
    print "=== E-VALUE DISTRIBUTION ==="
    printf "  Very significant (<1e-50):  %d (%.1f%%)\n", very_sig, 100*very_sig/total_hits
    printf "  Significant (1e-10–1e-50):  %d (%.1f%%)\n", sig, 100*sig/total_hits
    printf "  Weak (1e-5–1e-10):          %d (%.1f%%)\n", weak, 100*weak/total_hits
    printf "  Insignificant (>1e-5):      %d (%.1f%%)\n", insignificant, 100*insignificant/total_hits
    
    # Alignment length
    print ""
    print "=== ALIGNMENT LENGTH DISTRIBUTION ==="
    printf "Average alignment length:  %.1f aa\n", aln_len_sum/total_hits
    print ""
    printf "  Short (<50 aa):           %d (%.1f%%)\n", short_aln, 100*short_aln/total_hits
    printf "  Medium (50-200 aa):       %d (%.1f%%)\n", med_aln, 100*med_aln/total_hits
    printf "  Long (>200 aa):           %d (%.1f%%)\n", long_aln, 100*long_aln/total_hits
    
    # Paralog statistics
    print ""
    print "=== PARALOG ANALYSIS ==="
    
    # Calculate paralog statistics
    max_paralogs = 0
    max_gene = ""
    paralog_sum = 0
    paralog_total = 0
    
    for (gene in paralog_count) {
        paralog_total++
        paralog_sum += paralog_count[gene]
        
        if (paralog_count[gene] > max_paralogs) {
            max_paralogs = paralog_count[gene]
            max_gene = gene
        }
    }
    
    avg_paralogs = paralog_sum / paralog_total
    
    printf "Genes with paralogs:       %d (%.1f%% of queries)\n", 
           paralog_total, 100*paralog_total/length(queries)
    printf "Average paralogs per gene: %.1f\n", avg_paralogs
    printf "Max paralogs (gene):       %d (%s)\n", max_paralogs, max_gene
    
    # Distribution of paralog counts
    for (gene in paralog_count) {
        pc = paralog_count[gene]
        if (pc == 1) single++
        else if (pc <= 5) small_fam++
        else if (pc <= 20) med_fam++
        else if (pc <= 100) large_fam++
        else huge_fam++
    }
    
    print ""
    print "Gene family size distribution:"
    printf "  Singletons (1 paralog):    %d (%.1f%%)\n", single, 100*single/paralog_total
    printf "  Small (2-5):               %d (%.1f%%)\n", small_fam, 100*small_fam/paralog_total
    printf "  Medium (6-20):             %d (%.1f%%)\n", med_fam, 100*med_fam/paralog_total
    printf "  Large (21-100):            %d (%.1f%%)\n", large_fam, 100*large_fam/paralog_total
    printf "  Huge (>100):               %d (%.1f%%)\n", huge_fam, 100*huge_fam/paralog_total
    
    # Write paralog counts to file
    print "gene_id", "paralog_count" > paralog_file
    for (gene in paralog_count) {
        print gene, paralog_count[gene] > paralog_file
    }
    
    # Find top 50 gene families
    print ""
    print "=== TOP 50 GENE FAMILIES ==="
    print "(Genes with most paralogs, see " top_file " for full list)"
    print ""
    
    # Sort by paralog count (descending)
    n = asorti(paralog_count, sorted_genes, "@val_num_desc")
    
    print "Rank", "Gene_ID", "Paralog_Count" > top_file
    
    display_count = (n < 50) ? n : 50
    for (i = 1; i <= display_count; i++) {
        gene = sorted_genes[i]
        count = paralog_count[gene]
        printf "%3d  %-30s  %d\n", i, gene, count
        print i, gene, count > top_file
    }
    
    # Continue writing rest to file only
    for (i = display_count + 1; i <= n; i++) {
        gene = sorted_genes[i]
        count = paralog_count[gene]
        print i, gene, count > top_file
    }
    
    print ""
    print "========================================================"
    print ""
    print "Analysis complete!"
    print ""
    print "Output files:"
    print "  Summary:           " FILENAME ".summary.txt"
    print "  Paralog counts:    " paralog_file
    print "  Top gene families: " top_file
    print ""
}
' "$INPUT_FILE"

echo "Results saved to: $OUTPUT_DIR"

exit 0
