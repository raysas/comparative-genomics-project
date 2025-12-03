#!/bin/bash

# --------------------------------------------------------------------
# -- Create pairwise FASTA files from pre-generated pairs list
# -- Ultra-fast parallel processing
# --------------------------------------------------------------------

set -euo pipefail

# -- default parameters
PAIRS_FILE=''
PROTEIN_FASTA=''
OUTPUT_DIR=''
NUM_JOBS=$(nproc)
BATCH_SIZE=1000

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# -- arguments
while getopts "p:f:o:j:b:h" flag; do
    case "${flag}" in
        p) PAIRS_FILE="${OPTARG}" ;;
        f) PROTEIN_FASTA="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        j) NUM_JOBS="${OPTARG}" ;;
        b) BATCH_SIZE="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -p FILE       Input pairs list file (pairs_list.tsv)
  -f FILE       Protein FASTA file (for sequence lookup)
  -o DIR        Output directory for pairwise FASTA files
  -j NUM        Number of parallel jobs (default: all CPUs)
  -b NUM        Batch size for processing (default: 1000)
  -h            Show this help

INPUT FORMAT:
  The pairs_list.tsv should have columns:
    family_dir | gene1 | gene2 | seq1 | seq2
  OR
    family | gene1 | gene2 (sequences will be looked up from FASTA)

EXAMPLES:
  # With sequences in the TSV file
  $0 -p pairs_list.tsv -o output/families/

  # With sequence lookup from FASTA
  $0 -p pairs_list.tsv -f proteins.fa -o output/families/ -j 32

EOF
            exit 0
            ;;
        *)
            echo "Invalid option. Use -h for help."
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$PAIRS_FILE" ]; then
    echo -e "${RED}ERROR: Pairs file required (-p)${NC}"
    echo "Use -h for help"
    exit 1
fi

if [ -z "$OUTPUT_DIR" ]; then
    echo -e "${RED}ERROR: Output directory required (-o)${NC}"
    echo "Use -h for help"
    exit 1
fi

# Check if pairs file exists
if [ ! -f "$PAIRS_FILE" ]; then
    echo -e "${RED}ERROR: Pairs file not found: $PAIRS_FILE${NC}"
    exit 1
fi

# Convert to absolute paths
PAIRS_FILE=$(realpath "$PAIRS_FILE")
OUTPUT_DIR=$(realpath -m "$OUTPUT_DIR")
if [ -n "$PROTEIN_FASTA" ]; then
    PROTEIN_FASTA=$(realpath "$PROTEIN_FASTA")
fi

echo -e "${GREEN}===================================="
echo " PAIRWISE FASTA CREATION"
echo "====================================${NC}"
echo " Pairs file : $PAIRS_FILE"
if [ -n "$PROTEIN_FASTA" ]; then
    echo " Protein DB : $PROTEIN_FASTA"
fi
echo " Output dir : $OUTPUT_DIR"
echo " CPUs       : $NUM_JOBS"
echo " Batch size : $BATCH_SIZE"
echo -e "${GREEN}====================================${NC}"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check the format of the pairs file
echo -e "${YELLOW}Analyzing input file format...${NC}"
FIRST_LINE=$(head -n1 "$PAIRS_FILE")
NUM_COLS=$(echo "$FIRST_LINE" | awk -F'\t' '{print NF}')

if [ "$NUM_COLS" -ge 5 ]; then
    echo "  Format: TSV with sequences included (5+ columns)"
    HAS_SEQUENCES=true
elif [ "$NUM_COLS" -eq 3 ]; then
    echo "  Format: TSV with IDs only (3 columns)"
    HAS_SEQUENCES=false
    if [ -z "$PROTEIN_FASTA" ]; then
        echo -e "${RED}ERROR: Protein FASTA required (-f) when pairs file has no sequences${NC}"
        exit 1
    fi
else
    echo -e "${RED}ERROR: Invalid format. Expected 3 or 5+ columns, found $NUM_COLS${NC}"
    exit 1
fi

# Count total pairs
TOTAL_PAIRS=$(wc -l < "$PAIRS_FILE")
echo -e "${GREEN}✓ Found $TOTAL_PAIRS pairs to process${NC}"

# Create temp directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Process based on format
if [ "$HAS_SEQUENCES" = true ]; then
    echo -e "${YELLOW}Creating FASTA files from embedded sequences...${NC}"
    
    # Function to process pairs with sequences
    process_pair_with_seq() {
        local line="$1"
        local output_base="$2"
        
        # Parse TSV line
        IFS=$'\t' read -r family gene1 gene2 seq1 seq2 rest <<< "$line"
        
        # Clean family name
        family=$(echo "$family" | sed 's|.*/||' | tr -d '\r' | tr -d '[:space:]')
        
        # Create family directory
        local fam_dir="${output_base}/${family}"
        mkdir -p "$fam_dir"
        
        # Output file
        local out_file="${fam_dir}/${gene1}_${gene2}.fa"
        
        # Skip if exists
        if [ -f "$out_file" ] && [ -s "$out_file" ]; then
            return 0
        fi
        
        # Write FASTA
        {
            echo ">$gene1"
            echo "$seq1"
            echo ">$gene2"
            echo "$seq2"
        } > "$out_file"
    }
    
    export -f process_pair_with_seq
    export OUTPUT_DIR
    
    # Process in parallel
    cat "$PAIRS_FILE" | \
        parallel -j "$NUM_JOBS" \
                 --progress \
                 --pipe -N "$BATCH_SIZE" \
                 'while IFS= read -r line; do 
                     process_pair_with_seq "$line" "'"$OUTPUT_DIR"'"
                  done'
    
else
    echo -e "${YELLOW}Creating sequence index from FASTA...${NC}"
    
    # Build sequence index
    SEQ_INDEX="$TEMP_DIR/seq_index.txt"
    awk '/^>/{if(seq) print id "\t" seq; id=substr($0,2); seq=""} 
         !/^>/{seq=seq$0} 
         END{if(seq) print id "\t" seq}' "$PROTEIN_FASTA" > "$SEQ_INDEX"
    
    SEQ_COUNT=$(wc -l < "$SEQ_INDEX")
    echo -e "${GREEN}✓ Indexed $SEQ_COUNT sequences${NC}"
    
    echo -e "${YELLOW}Creating FASTA files with sequence lookup...${NC}"
    
    # Function to process pairs with lookup
    process_pair_with_lookup() {
        local family="$1"
        local gene1="$2"
        local gene2="$3"
        local seq_index="$4"
        local output_base="$5"
        
        # Clean family name
        family=$(echo "$family" | sed 's|family||' | tr -d '\r' | tr -d '[:space:]')
        
        # Create family directory
        local fam_dir="${output_base}/family${family}"
        mkdir -p "$fam_dir"
        
        # Output file
        local out_file="${fam_dir}/${gene1}_${gene2}.fa"
        
        # Skip if exists
        if [ -f "$out_file" ] && [ -s "$out_file" ]; then
            return 0
        fi
        
        # Look up sequences
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
    
    export -f process_pair_with_lookup
    export OUTPUT_DIR SEQ_INDEX
    
    # Process in parallel
    cat "$PAIRS_FILE" | \
        parallel -j "$NUM_JOBS" \
                 --progress --bar \
                 --colsep '\t' \
                 "process_pair_with_lookup {1} {2} {3} '$SEQ_INDEX' '$OUTPUT_DIR'"
fi

# Count results
echo ""
echo -e "${YELLOW}Counting created files...${NC}"
CREATED_FILES=$(find "$OUTPUT_DIR" -name "*.fa" -type f | wc -l)
CREATED_FAMILIES=$(find "$OUTPUT_DIR" -maxdepth 1 -type d -name "family*" | wc -l)

# Final report
echo ""
echo -e "${GREEN}===================================="
echo " COMPLETED"
echo "====================================${NC}"
echo " Total pairs    : $TOTAL_PAIRS"
echo " Files created  : $CREATED_FILES"
echo " Family folders : $CREATED_FAMILIES"
echo " Output location: $OUTPUT_DIR"
echo -e "${GREEN}====================================${NC}"

if [ "$CREATED_FILES" -lt "$TOTAL_PAIRS" ]; then
    MISSING=$((TOTAL_PAIRS - CREATED_FILES))
    echo -e "${YELLOW}Note: $MISSING pairs may have been skipped (already exist or missing sequences)${NC}"
fi

exit 0