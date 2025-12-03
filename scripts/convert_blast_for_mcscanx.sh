#!/bin/bash
# Convert BLAST tabular output to MCScanX format (4 columns: gene1 gene2 evalue score)
# Usage: bash scripts/convert_blast_to_mcscanx.sh input.tsv output.blast [species_prefix]

set -euo pipefail

INPUT="${1:-}"
OUTPUT="${2:-}"
PREFIX="${3:-}"

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Usage: $0 <input_blast.tsv> <output.blast> [species_prefix]"
    echo "  Example: $0 output/blast_output/blast_results_with_coverage.tsv output/mcscanx/soybean.blast soybean"
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT" >&2
    exit 1
fi

# Detect columns (handle both with/without header)
# Expected columns: qseqid, sseqid (or gene1, gene2), evalue, bitscore
# Try to auto-detect column positions

# Read first line to check for header
read -r FIRST_LINE < "$INPUT"

# Check if it looks like a header (contains common BLAST column names)
if [[ "$FIRST_LINE" =~ (qseqid|gene1|query) ]]; then
    HAS_HEADER=1
    # Parse header to find column indices
    IFS=$'\t' read -r -a COLS <<< "$FIRST_LINE"
    
    Q_COL=-1
    S_COL=-1
    E_COL=-1
    B_COL=-1
    
    for i in "${!COLS[@]}"; do
        col="${COLS[$i],,}"  # lowercase
        case "$col" in
            qseqid|query|qacc|gene1|gene) Q_COL=$((i+1)) ;;
            sseqid|subject|sacc|gene2) S_COL=$((i+1)) ;;
            evalue|e-val|expect) E_COL=$((i+1)) ;;
            bitscore|score) B_COL=$((i+1)) ;;
        esac
    done
    
    if [[ $Q_COL -eq -1 || $S_COL -eq -1 ]]; then
        echo "ERROR: Cannot find query/subject columns in header" >&2
        exit 1
    fi
    
    # Default to small values if evalue/score not found
    [[ $E_COL -eq -1 ]] && E_COL="1e-5"
    [[ $B_COL -eq -1 ]] && B_COL="50"
    
    SKIP_LINES=1
else
    # No header - assume standard BLAST -outfmt 6 order
    Q_COL=1
    S_COL=2
    E_COL=11
    B_COL=12
    SKIP_LINES=0
fi

echo "Converting BLAST to MCScanX format..."
echo "  Input: $INPUT"
echo "  Output: $OUTPUT"
[[ -n "$PREFIX" ]] && echo "  Prefix: $PREFIX"

# Process BLAST file
{
    if [[ $SKIP_LINES -eq 1 ]]; then
        tail -n +2 "$INPUT"  # Skip header
    else
        cat "$INPUT"
    fi
} | awk -F'\t' -v q="$Q_COL" -v s="$S_COL" -v e="$E_COL" -v b="$B_COL" -v pfx="$PREFIX" '
BEGIN {
    OFS="\t"
    # Handle constant values for missing columns
    if (e !~ /^[0-9]+$/) { e_const = e; e = -1 } else { e_const = "" }
    if (b !~ /^[0-9]+$/) { b_const = b; b = -1 } else { b_const = "" }
}
{
    # Skip self-hits
    if ($q == $s) next
    
    gene1 = $q
    gene2 = $s
    
    # Add prefix if provided
    if (pfx != "") {
        gene1 = pfx "|" gene1
        gene2 = pfx "|" gene2
    }
    
    # Get evalue and score
    if (e > 0) {
        evalue = $e
    } else {
        evalue = e_const
    }
    
    if (b > 0) {
        score = $b
    } else {
        score = b_const
    }
    
    # MCScanX format: gene1 gene2 evalue score
    printf "%s\t%s\t%s\t%s\n", gene1, gene2, evalue, score
}' > "$OUTPUT"

# Count lines
LINE_COUNT=$(wc -l < "$OUTPUT")
echo "✓ Converted $LINE_COUNT BLAST pairs to MCScanX format"
echo "  Saved: $OUTPUT"
