#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    Calculates Ks/Ka values using PAML yn00 from codon alignments
#    (processes back-translated alignments from step 4)
# -- Usage:
#    bash ./pipeline_2/5_calculate_ks.sh [-i INPUT_DIR] [-o OUTPUT_DIR] [-h]
# -- default (without params) equivalent to:
#    Auto-detects from species directories
# --------------------------------------------------------------------

# -- default parameters
INPUT_DIR=''
OUTPUT_DIR=''
AUTO_DETECT=true

# -- arguments
while getopts "i:o:h" flag; do
    case "${flag}" in
        i) INPUT_DIR="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i DIR        Input directory with codon alignments
  -o DIR        Output directory for Ks/Ka results
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
    Input:  output/{species}/codon_alignments/
    Output: output/{species}/ks_results/

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i output/glycine_max/codon_alignments/ \\
     -o output/glycine_max/ks_results/

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
if [ -n "$INPUT_DIR" ] || [ -n "$OUTPUT_DIR" ]; then
    AUTO_DETECT=false
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    SPECIES_DIRS=(output/*/codon_alignments)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        INPUT_DIR="${SPECIES_DIRS[0]}"
        SPECIES_NAME=$(basename $(dirname "$INPUT_DIR"))
        OUTPUT_DIR="output/${SPECIES_NAME}/ks_results"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Input:  $INPUT_DIR"
        echo "   Output: $OUTPUT_DIR"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species codon alignment directories found"
        echo "       Specify paths explicitly with -i and -o"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No codon alignment results found in output/*/"
        echo "       Run pipeline_2/4_backtranslate.sh first"
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

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    echo "       Run pipeline_2/4_backtranslate.sh first"
    exit 1
fi

# Check if yn00 is available
if ! command -v yn00 &> /dev/null; then
    echo "ERROR: yn00 (PAML) not found. Please install:"
    echo "       conda install -c bioconda paml"
    echo "   or  download from: http://abacus.gene.ucl.ac.uk/software/paml.html"
    exit 1
fi

# Check for yn00 control template
YN00_TEMPLATE="$(dirname "${BASH_SOURCE[0]}")/yn00.ctl_master"
if [ ! -f "$YN00_TEMPLATE" ]; then
    echo "Creating default yn00 control template..."
    mkdir -p "$(dirname "$YN00_TEMPLATE")"
    cat > "$YN00_TEMPLATE" <<EOF
seqfile = SEQFILE_PLACEHOLDER * sequence data file name
outfile = OUTFILE_PLACEHOLDER * main result file
verbose = 0  * 1: detailed output (list sequences), 0: concise output
icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
weighting = 0  * weighting pathways between codons (0/1)?
commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)?
ndata = 1
 * Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
 * 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,
* 10: blepharisma nu.
 * These codes correspond to transl_table 1 to 11 of GENEBANK.
EOF
fi

echo "-- Calculating Ks/Ka values with PAML yn00..."

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Count codon alignments for progress
TOTAL_ALIGNMENTS=$(find "$INPUT_DIR" -name "*_codon.aln" | wc -l)

echo "   Found $TOTAL_ALIGNMENTS codon alignments for Ks calculation"

if [ $TOTAL_ALIGNMENTS -eq 0 ]; then
    echo "ERROR: No codon alignment files (*_codon.aln) found in $INPUT_DIR"
    echo "       Run pipeline_2/4_backtranslate.sh first"
    exit 1
fi

START_TIME=$(date +%s)
PROCESSED=0
SUCCESS=0
FAILED=0

# Process all family directories
for fam_dir in "$INPUT_DIR"/family*; do
    if [ ! -d "$fam_dir" ]; then
        continue
    fi
    
    fam_name=$(basename "$fam_dir")
    fam_out="$OUTPUT_DIR/$fam_name"
    mkdir -p "$fam_out"

    echo "▶ Processing $fam_name"

    # Process all codon alignment files in this family
    for aln_file in "$fam_dir"/*_codon.aln; do
        if [ ! -f "$aln_file" ]; then
            continue
        fi
        
        aln_base=$(basename "$aln_file" _codon.aln)
        
        # Create control file for this alignment
        ctl_file="${fam_out}/${aln_base}_yn00.ctl"
        out_file="${fam_out}/${aln_base}_yn00.out"
        
        # Create yn00 control file from template (use absolute paths)
        abs_aln_file=$(realpath "$aln_file")
        abs_out_file=$(cd "$(dirname "$out_file")" && pwd)/$(basename "$out_file")
        sed -e "s|SEQFILE_PLACEHOLDER|$abs_aln_file|g" \
            -e "s|OUTFILE_PLACEHOLDER|$abs_out_file|g" \
            "$YN00_TEMPLATE" > "$ctl_file"
        
        # Run yn00
        # Run yn00
        if (cd "$(dirname "$ctl_file")" && rm -f yn00.ctl && yn00 "$(basename "$ctl_file")") >/dev/null 2>&1; then
            # Check if output file was created and contains results
            if [ -f "$out_file" ] && grep -q "dN/dS" "$out_file" 2>/dev/null; then
                echo "   ✓ $aln_base -> ${aln_base}_yn00.out"
                ((SUCCESS++))
            else
                echo "   ✗ $aln_base: yn00 produced no results"
                rm -f "$out_file"
                ((FAILED++))
            fi
        else
            echo "   ✗ $aln_base: yn00 failed"
            rm -f "$out_file"
            ((FAILED++))
        fi
        
        # Clean up temporary files (yn00 creates these in the working directory)
        (cd "$(dirname "$ctl_file")" && rm -f rst rst1 2nstrees 4fold.nuc lnf yn00.ctl)
        rm -f "$ctl_file"
        
        ((PROCESSED++))
        
        # Progress update every 50 alignments
        if [ $((PROCESSED % 50)) -eq 0 ]; then
            ELAPSED=$(($(date +%s) - START_TIME))
            echo "   Progress: $PROCESSED/$TOTAL_ALIGNMENTS alignments (${ELAPSED}s elapsed)"
        fi
    done
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "-- Ks/Ka calculation completed"
echo "   Total processed: $PROCESSED"
echo "   Successful: $SUCCESS"
echo "   Failed: $FAILED"
echo "   Processing time: ${ELAPSED}s"
echo "   Output directory: $OUTPUT_DIR"

if [ $SUCCESS -gt 0 ]; then
    echo ""
    echo "-- Next step: Consolidate results"
    echo "   bash pipeline_2/6_consolidate_results.sh"
fi

exit 0