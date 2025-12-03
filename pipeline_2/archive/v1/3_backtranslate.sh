#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    Back-translates protein alignments to codon alignments using pal2nal
#    (converts ClustalW .aln files to codon-based alignments)
# -- Usage:
#    bash ./pipeline_2/4_backtranslate.sh [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CDS_FILE] [-h]
# -- default (without params) equivalent to:
#    Auto-detects from species directories
# --------------------------------------------------------------------

# -- default parameters
INPUT_DIR=''
OUTPUT_DIR=''
CDS_FILE=''
AUTO_DETECT=true

# -- arguments
while getopts "i:o:c:h" flag; do
    case "${flag}" in
        i) INPUT_DIR="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        c) CDS_FILE="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i DIR        Input directory with protein alignments (.aln files)
  -o DIR        Output directory for codon alignments
  -c FILE       CDS FASTA file for backtranslation
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
    Input:  output/{species}/alignments/
    Output: output/{species}/codon_alignments/
    CDS:    data/{species}/cds.fa

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i output/glycine_max/alignments/ \\
     -o output/glycine_max/codon_alignments/ \\
     -c data/glycine_max/cds.fa

EOF
            exit 0
            ;;
        *)
            echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
    esac
done

# If explicit arguments provided, disable auto-detection
if [ -n "$INPUT_DIR" ] || [ -n "$OUTPUT_DIR" ] || [ -n "$CDS_FILE" ]; then
    AUTO_DETECT=false
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    SPECIES_DIRS=(output/*/alignments)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        INPUT_DIR="${SPECIES_DIRS[0]}"
        SPECIES_NAME=$(basename $(dirname "$INPUT_DIR"))
        OUTPUT_DIR="output/${SPECIES_NAME}/codon_alignments"
        CDS_FILE="data/${SPECIES_NAME}/cds.fa"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Input:  $INPUT_DIR"
        echo "   Output: $OUTPUT_DIR"
        echo "   CDS:    $CDS_FILE"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species alignment directories found"
        echo "       Specify paths explicitly with -i, -o, -c"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No alignment results found in output/*/"
        echo "       Run pipeline_2/3_align_proteins.sh first"
        exit 1
    fi
fi

# Setup logging
LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo "-- Parameters:"
echo "   INPUT DIR  : $INPUT_DIR"
echo "   OUTPUT DIR : $OUTPUT_DIR"
echo "   CDS FILE   : $CDS_FILE"

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    echo "       Run pipeline_2/3_align_proteins.sh first"
    exit 1
fi

# Check CDS file exists
if [ ! -f "$CDS_FILE" ]; then
    echo "ERROR: CDS file not found: $CDS_FILE"
    echo "       Run pipeline_1/0_extract_data.sh first"
    exit 1
fi

# Check if pal2nal.pl exists (try multiple locations)
PAL2NAL_SCRIPT=""
if [ -f "pipeline_2/pal2nal.pl" ]; then
    PAL2NAL_SCRIPT="pipeline_2/pal2nal.pl"
elif [ -f "pal2nal.pl" ]; then
    PAL2NAL_SCRIPT="pal2nal.pl"
elif [ -f "$(dirname "$0")/pal2nal.pl" ]; then
    PAL2NAL_SCRIPT="$(dirname "$0")/pal2nal.pl"
fi

if [ -z "$PAL2NAL_SCRIPT" ] || [ ! -f "$PAL2NAL_SCRIPT" ]; then
    echo "ERROR: pal2nal.pl not found. Checked locations:"
    echo "       pipeline_2/pal2nal.pl"
    echo "       pal2nal.pl"
    echo "       $(dirname "$0")/pal2nal.pl"
    echo "       Download from: http://www.bork.embl.de/pal2nal/"
    exit 1
fi

echo "-- Back-translating protein alignments to codon alignments..."

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count alignments for progress
TOTAL_ALIGNMENTS=$(find "$INPUT_DIR" -name "*.aln" | wc -l)

echo "   Found $TOTAL_ALIGNMENTS protein alignments to back-translate"

if [ $TOTAL_ALIGNMENTS -eq 0 ]; then
    echo "ERROR: No alignment files (.aln) found in $INPUT_DIR"
    echo "       Run pipeline_2/3_align_proteins.sh first"
    exit 1
fi

START_TIME=$(date +%s)
PROCESSED=0
SUCCESS=0
FAILED=0

# Load CDS sequences for fast lookup
echo "   Loading CDS sequences..."
awk '/^>/{if(seq) printf "%s\t%s\n", seq_id, seq; seq_id=substr($0,2); seq=""} !/^>/{seq=seq$0} END{if(seq) printf "%s\t%s\n", seq_id, seq}' "$CDS_FILE" > "${OUTPUT_DIR}/.cds_lookup.tmp"

# Process all family directories
for fam_dir in "$INPUT_DIR"/family*; do
    if [ ! -d "$fam_dir" ]; then
        continue
    fi
    
    fam_name=$(basename "$fam_dir")
    fam_out="$OUTPUT_DIR/$fam_name"
    mkdir -p "$fam_out"

    echo "▶ Processing $fam_name"

    # Process all alignment files in this family
    for aln_file in "$fam_dir"/*.aln; do
        if [ ! -f "$aln_file" ]; then
            continue
        fi
        
        aln_base=$(basename "$aln_file" .aln)
        
        # Extract gene IDs from filename (gene1_gene2.aln)
        gene1=$(echo "$aln_base" | cut -d'_' -f1)
        gene2=$(echo "$aln_base" | cut -d'_' -f2)
        
        # Create temporary CDS file for this pair
        temp_cds="${fam_out}/${aln_base}_temp.fa"
        
        # Extract CDS sequences for this pair
        grep -E "^${gene1}|^${gene2}" "${OUTPUT_DIR}/.cds_lookup.tmp" | \
        awk -F'\t' '{printf ">%s\n%s\n", $1, $2}' > "$temp_cds"
        
        # Check if both sequences were found
        cds_count=$(grep -c "^>" "$temp_cds")
        if [ $cds_count -ne 2 ]; then
            echo "   ✗ $aln_base: Missing CDS sequences ($cds_count/2 found)"
            rm -f "$temp_cds"
            ((FAILED++))
            ((PROCESSED++))
            continue
        fi
        
        # Back-translate with pal2nal
        output_file="${fam_out}/${aln_base}_codon.aln"
        
        if perl "$PAL2NAL_SCRIPT" "$aln_file" "$temp_cds" -output paml -nogap > "$output_file" 2>/dev/null; then
            # echo "   ✓ $aln_base -> ${aln_base}_codon.aln"
            ((SUCCESS++))
        else
            # echo "   ✗ $aln_base: pal2nal failed"
            rm -f "$output_file"
            ((FAILED++))
        fi
        
        # Cleanup temp file
        rm -f "$temp_cds"
        ((PROCESSED++))
        
        # Progress update every 100 alignments
        if [ $((PROCESSED % 100)) -eq 0 ]; then
            ELAPSED=$(($(date +%s) - START_TIME))
            echo "   Progress: $PROCESSED/$TOTAL_ALIGNMENTS alignments (${ELAPSED}s elapsed)"
        fi
    done
done

# Cleanup
rm -f "${OUTPUT_DIR}/.cds_lookup.tmp"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "-- Back-translation completed"
echo "   Total processed: $PROCESSED"
echo "   Successful: $SUCCESS"
echo "   Failed: $FAILED"
echo "   Processing time: ${ELAPSED}s"
echo "   Output directory: $OUTPUT_DIR"

if [ $SUCCESS -gt 0 ]; then
    echo ""
    echo "-- Next step: Calculate Ks values"
    echo "   bash pipeline_2/5_calculate_ks.sh"
fi

exit 0