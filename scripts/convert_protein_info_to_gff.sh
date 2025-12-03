#!/bin/bash
# Convert protein_info_longest.csv to MCScanX GFF format
# Format: gene_id chromosome start end
# Usage: bash scripts/convert_protein_info_to_gff.sh [input.csv] [output.gff] [species_prefix]

set -euo pipefail

INPUT="${1:-data/protein_info_longest.csv}"
OUTPUT="${2:-output/mcscanx/soybean.gff}"
PREFIX="${3:-gm}"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input file not found: $INPUT" >&2
    exit 1
fi

echo "Converting protein info to MCScanX GFF format..."
echo "  Input: $INPUT"
echo "  Output: $OUTPUT"
echo "  Chromosome prefix: $PREFIX"

# Create output directory if needed
mkdir -p "$(dirname "$OUTPUT")"

# Convert CSV to GFF format
# Skip header, extract: peptide_id (matches BLAST), chromosome, start_pos, end_pos
# Add chromosome prefix (e.g., "1" -> "gm1")
# Handle both integer and float chromosome values (1.0 -> 1)

awk -F',' 'NR>1 {
    gene_id = $1        # peptide_id (e.g., KRH*, KRG*)
    chr = $5            # chromosome
    start = $6          # start_pos
    end = $7            # end_pos
    
    # Skip if any field is empty or "nan"
    if (gene_id == "" || chr == "" || chr == "nan" || start == "" || end == "") next
    
    # Convert float chromosome to int (1.0 -> 1)
    chr = int(chr)
    if (chr < 1) next
    
    # Convert start/end to int
    start = int(start)
    end = int(end)
    if (start <= 0 || end <= 0) next
    
    # Format: chr_with_prefix gene_id start end
    printf "%s\t%s\t%d\t%d\n", "'$PREFIX'"chr, gene_id, start, end
}' "$INPUT" > "$OUTPUT"

LINE_COUNT=$(wc -l < "$OUTPUT")
echo "  Converted $LINE_COUNT genes"

# Show sample
# echo ""
# echo "Sample output (first 5 lines):"
# head -5 "$OUTPUT"
# echo ""
# echo "✓ GFF file created: $OUTPUT"
