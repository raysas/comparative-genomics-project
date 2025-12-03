#!/bin/bash

# handles command line length limitation

# -- default parameters
FAMILY_FILE=''
PROTEIN_FASTA=''
OUTPUT_DIR=''
AUTO_DETECT=true
NUM_JOBS=$(nproc)

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "-- Loading protein sequences into indexed file..."
# Create an indexed sequence file for fast lookup
SEQ_INDEX="/tmp/protein_sequences_$$.idx"
awk '/^>/{if(seq) printf "%s\t%s\n", id, seq; id=substr($0,2); seq=""} 
     !/^>/{seq=seq$0} 
     END{if(seq) printf "%s\t%s\n", id, seq}' "$PROTEIN_FASTA" > "$SEQ_INDEX"

echo "-- Generating pair list (IDs only, no sequences)..."
PAIR_LIST="/tmp/pair_list_$$.txt"

# Generate list of pairs (just IDs, not sequences)
awk -F'\t' 'NR > 1 {
    family = $2
    gene = $1
    gsub(/[[:space:]\r]/, "", family)
    gsub(/[[:space:]\r]/, "", gene)
    genes[family] = genes[family] ? genes[family] "," gene : gene
}
END {
    for (fam in genes) {
        split(genes[fam], gene_list, ",")
        n = length(gene_list)
        for (i = 1; i <= n; i++) {
            for (j = i+1; j <= n; j++) {
                print "family" fam, gene_list[i], gene_list[j]
            }
        }
    }
}' "$FAMILY_FILE" > "$PAIR_LIST"

TOTAL_PAIRS=$(wc -l < "$PAIR_LIST")
echo "-- Creating $TOTAL_PAIRS pairwise FASTA files..."

# Function that looks up sequences from the index file
make_fasta_from_index() {
    local fam_dir="$1"
    local gene1="$2"
    local gene2="$3"
    local seq_index="$4"
    local output_base="$5"
    
    # Create family directory if needed
    mkdir -p "$output_base/$fam_dir"
    
    local out_file="$output_base/$fam_dir/${gene1}_${gene2}.fa"
    
    # Skip if exists
    if [ -f "$out_file" ] && [ -s "$out_file" ]; then
        return 0
    fi
    
    # Look up sequences from index
    seq1=$(grep "^${gene1}	" "$seq_index" | cut -f2)
    seq2=$(grep "^${gene2}	" "$seq_index" | cut -f2)
    
    if [ -n "$seq1" ] && [ -n "$seq2" ]; then
        {
            echo ">$gene1"
            echo "$seq1"
            echo ">$gene2"
            echo "$seq2"
        } > "$out_file"
    fi
}
export -f make_fasta_from_index

# Process with parallel (no sequences in command line)
cat "$PAIR_LIST" | \
    parallel -j "$NUM_JOBS" --progress --colsep ' ' \
    "make_fasta_from_index {1} {2} {3} '$SEQ_INDEX' '$OUTPUT_DIR'"

# Cleanup
rm -f "$SEQ_INDEX" "$PAIR_LIST"

echo "-- Completed creating pairwise FASTA files"