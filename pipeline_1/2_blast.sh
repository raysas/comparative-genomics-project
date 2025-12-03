#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    1) creates a BLAST database from a given peptide fasta file
#    2) runs a BLASTP search of the peptides against themselves
# !! if the output directory already exists, it will skip running blast unless -f flag is provided
# -- Usage:
#    bash ./pipeline_1/2_blast.sh [-i INPUT_FILE] [-o OUTPUT_DIR] [-f] [-h]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/2_blast.sh -i "data/peptides_longest.fa" -o "output/blast_output"
# --------------------------------------------------------------------
###########################################################################

# ---------------------------------------------------------------------
# prepare variables, get arguments, set up logging
# ---------------------------------------------------------------------
# -- message on what this script does
# cat <<EOF
# -- this script:
#   1) creates a BLAST database from a given peptide fasta file
#   2) runs a BLASTP search of the peptides against themselves
# EOF

# Resolve paths to run from anywhere and improve UX
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
NC='\033[0m'


# -- default parameters
INPUT_FILE=''
OUTPUT_DIR=''
AUTO_DETECT=true
E_VALUE='10'
MAX_TARGET_SEQS=500

DB_NAME='peptide_db'
OUTPUT_FILE='blast_results.tsv'

# -- make a flag to not keep blast_output directory across runs and make a new one each time
force_new=false

# -- get arguments
SPECIES_NAME=""
while getopts "i:o:s:e:fh" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}";;
        o) OUTPUT_DIR="${OPTARG}";;
        s) SPECIES_NAME="${OPTARG}"; AUTO_DETECT=true;;
        e) E_VALUE="${OPTARG}";;
        f) force_new=true ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE       Input FASTA file (peptides_longest.fa)
  -o DIR        Output directory for BLAST results
    -s NAME       Species name (uses data/NAME and output/NAME paths)
    -e VALUE      E-value threshold for BLAST (default: 1e-5)
  -f            Force new output directory (rename existing)
  -h            Show this help

AUTO-DETECTION:
    If -s is provided or no options specified, automatically detects species directory:
        Input :  data/{species}/peptides_longest.fa
        Output: output/pipeline1/{species}/blast_results/

EXAMPLES:
  # Auto-detect (recommended)
  $0

    # Explicit paths
    $0 -i data/glycine_max/peptides_longest.fa -o output/glycine_max/blast_output

    # Species-driven auto paths
    $0 -s glycine_max

  # Force new run (archive existing results)
  $0 -f

EOF
            exit 0
            ;;
        *)
            echo "Invalid option: -${OPTARG}" >&2
            exit 1
            ;;
    esac
done

# Check if user provided partial input - disable auto-detect only if both are specified
if [ -n "$INPUT_FILE" ] || [ -n "$OUTPUT_DIR" ]; then
    AUTO_DETECT=false
    
    # If only one is specified, auto-detect the other
    if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
        echo "-- Partial arguments detected, auto-detecting missing paths..."
        
        SPECIES_DIRS=(data/*/)
        if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
            SPECIES_DIR="${SPECIES_DIRS[0]%/}"
            SPECIES_NAME=$(basename "$SPECIES_DIR")
            
            [ -z "$INPUT_FILE" ] && INPUT_FILE="${SPECIES_DIR}/processed/peptides_longest.fa"
            [ -z "$OUTPUT_DIR" ] && OUTPUT_DIR="output/pipeline1/${SPECIES_NAME}/blast_results"
            
            echo "   Detected species: $SPECIES_NAME"
        else
            echo "ERROR: Cannot auto-detect. Provide both -i and -o"
            exit 1
        fi
    fi
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo -e "${YELLOW}-- Auto-detecting species directory...${NC}"
    if [ -n "$SPECIES_NAME" ]; then
        SPECIES_DIR="$REPO_ROOT/data/$SPECIES_NAME"
        if [ -d "$SPECIES_DIR" ]; then
            INPUT_FILE="${SPECIES_DIR}/processed/peptides_longest.fa"
            OUTPUT_DIR="${REPO_ROOT}/output/pipeline1/${SPECIES_NAME}/blast_results"
            echo -e "   Detected species: ${GREEN}$SPECIES_NAME${NC}"
            echo -e "   Input : ${BLUE}$INPUT_FILE${NC}"
            echo -e "   Output: ${BLUE}$OUTPUT_DIR${NC}"
        else
            echo -e "${RED}ERROR:${NC} Species directory not found: $SPECIES_DIR"
            exit 1
        fi
    else
        SPECIES_DIRS=("$REPO_ROOT"/data/*/)
        if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
            SPECIES_DIR="${SPECIES_DIRS[0]%/}"
            SPECIES_NAME=$(basename "$SPECIES_DIR")
            INPUT_FILE="${SPECIES_DIR}/processed/peptides_longest.fa"
            OUTPUT_DIR="${REPO_ROOT}/output/pipeline1/${SPECIES_NAME}/blast_results"
            echo -e "   Detected species: ${GREEN}$SPECIES_NAME${NC}"
            echo -e "   Input : ${BLUE}$INPUT_FILE${NC}"
            echo -e "   Output: ${BLUE}$OUTPUT_DIR${NC}"
        elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
            echo -e "${RED}ERROR: Multiple species directories found in data/. Use -s, -i and -o.${NC}"
            exit 1
        else
            echo -e "${RED}ERROR: No species directories found and no default inputs present.${NC}"
            echo "       Run 1_filter_isoforms.sh first or provide -i and -o."
            exit 1
        fi
    fi
fi


# -- name log file based on script name
# Logging inside pipeline directory
LOG_DIR="${SCRIPT_DIR}/logs/pipeline"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh)_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo -e "${GREEN}====================================${NC}"
echo -e "${GREEN} Step 2: BLAST ${NC}"
echo -e "${GREEN}====================================${NC}"
echo -e "Parameters:"
echo -e "  INPUT  : ${BLUE}$INPUT_FILE${NC}"
echo -e "  OUTPUT : ${BLUE}$OUTPUT_DIR${NC}"
echo -e "  E-VALUE: ${YELLOW}$E_VALUE${NC}"
echo -e "  FORCE  : ${YELLOW}$force_new${NC}"

# -- check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    echo "       Run 1_filter_isoforms.sh first"
    exit 1
fi

# -------------------------------------------------------------------------
# -- main script logic
# -------------------------------------------------------------------------

# -- handle existing output directory
if [ "$force_new" = false ] && [ -d "$OUTPUT_DIR" ]; then
    echo "-- output directory $OUTPUT_DIR already exists and -f flag not set; using existing directory"
    echo "-- to force creation of a new output directory, use the -f flag"
# -- else output directory exists and -f flag is set: rename exsisting output dire to old_blast_output_1, old_blast_output_2, etc.
else
    echo -e "${YELLOW}-- Preparing output directory:${NC} $OUTPUT_DIR"
    leaf_dir=$(basename "$OUTPUT_DIR")
    parent_dir=$(dirname "$OUTPUT_DIR")
    if [ -d "$OUTPUT_DIR" ]; then
        # -- move leaf to leaf_old_numbered
        n=1
        while [ -d "${parent_dir}/${leaf_dir}_old_${n}" ]; do
            ((n++))
        done
        mv "$OUTPUT_DIR" "${parent_dir}/${leaf_dir}_old_${n}"
        echo -e "${YELLOW}-- Moved existing:${NC} $OUTPUT_DIR -> ${parent_dir}/${leaf_dir}_old_${n}"
    fi

    # -- create output directory
    mkdir -p "$OUTPUT_DIR"
    echo -e "${GREEN}✓ Created output directory:${NC} $OUTPUT_DIR"

    # -- make blast db
    if ! command -v makeblastdb &> /dev/null
    then
        echo -e "${RED}ERROR:${NC} makeblastdb not found. Install BLAST+ tools."
        exit 2
    fi

    makeblastdb -in "$INPUT_FILE" -dbtype prot -out "${OUTPUT_DIR}/${DB_NAME}"
    echo -e "${GREEN}✓ BLAST DB created from:${NC} $INPUT_FILE"
    
    # Count total sequences for progress estimation
    TOTAL_SEQS=$(grep -c "^>" "$INPUT_FILE")
    echo -e "Total sequences to process: ${YELLOW}$TOTAL_SEQS${NC}"

    # -- run blastp with optimizations for speed
    if ! command -v blastp &> /dev/null
    then
        echo -e "${RED}ERROR:${NC} blastp not found. Install BLAST+ tools."
        exit 3
    fi
    
    # Detect available CPU cores
    NCORES=$(nproc 2>/dev/null || echo 4)
    # Use all cores - 1
    BLAST_THREADS=$((NCORES - 1))
    [ $BLAST_THREADS -lt 1 ] && BLAST_THREADS=1
    
    echo -e "${YELLOW}Running BLASTP with${NC} $BLAST_THREADS threads (detected $NCORES cores)"
    echo "   Progress will be logged to: ${OUTPUT_DIR}/blast_progress.txt"
    echo ""
    
    # Create a progress monitoring script
    PROGRESS_LOG="${OUTPUT_DIR}/blast_progress.txt"
    TEMP_OUTPUT="${OUTPUT_DIR}/.blast_results_temp.tsv"
    
    # Start timestamp
    START_TIME=$(date +%s)
    echo "Started: $(date)" > "$PROGRESS_LOG"
    
    # Run BLAST in background and monitor progress
        blastp -query "$INPUT_FILE" -db "${OUTPUT_DIR}/${DB_NAME}" \
            -out "$TEMP_OUTPUT" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -evalue "$E_VALUE" \
            -max_target_seqs $MAX_TARGET_SEQS \
            -num_threads $BLAST_THREADS \
            -word_size 6 \
            -threshold 21 \
            -comp_based_stats 0 \
            -seg no &
    
    BLAST_PID=$!
    
    echo -e "BLAST process started (PID: ${YELLOW}$BLAST_PID${NC})"
    echo ""
    echo "Monitoring progress (press Ctrl+C to stop monitoring, BLAST will continue):"
    echo "───────────────────────────────────────────────────────────────────"
    
    # Monitor progress while BLAST runs
    while kill -0 $BLAST_PID 2>/dev/null; do
        if [ -f "$TEMP_OUTPUT" ]; then
            CURRENT_HITS=$(wc -l < "$TEMP_OUTPUT" 2>/dev/null || echo 0)
            ELAPSED=$(($(date +%s) - START_TIME))
            
            # Estimate progress (rough approximation)
            QUERIES_DONE=$(cut -f1 "$TEMP_OUTPUT" 2>/dev/null | sort -u | wc -l || echo 0)
            
            if [ $TOTAL_SEQS -gt 0 ] && [ $QUERIES_DONE -gt 0 ]; then
                PERCENT=$((QUERIES_DONE * 100 / TOTAL_SEQS))
                [ $PERCENT -gt 100 ] && PERCENT=100
                
                # Estimate time remaining
                if [ $PERCENT -gt 0 ]; then
                    TOTAL_EST=$((ELAPSED * 100 / PERCENT))
                    REMAINING=$((TOTAL_EST - ELAPSED))
                    
                    printf "\r  Progress: %3d%% | Queries: %d/%d | Hits: %s | Elapsed: %02d:%02d:%02d | ETA: %02d:%02d:%02d" \
                        $PERCENT $QUERIES_DONE $TOTAL_SEQS $CURRENT_HITS \
                        $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60)) \
                        $((REMAINING/3600)) $(((REMAINING%3600)/60)) $((REMAINING%60))
                else
                    printf "\r  Processing... | Hits: %s | Elapsed: %02d:%02d:%02d" \
                        $CURRENT_HITS \
                        $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
                fi
            else
                printf "\r  Starting... | Hits: %s | Elapsed: %02d:%02d:%02d" \
                    $CURRENT_HITS \
                    $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
            fi
            
            # Log progress to file
            echo "$(date '+%Y-%m-%d %H:%M:%S') | Queries: $QUERIES_DONE/$TOTAL_SEQS | Hits: $CURRENT_HITS" >> "$PROGRESS_LOG"
        fi
        sleep 5
    done
    
    # Wait for BLAST to finish
    wait $BLAST_PID
    BLAST_EXIT=$?
    
    echo ""
    echo "───────────────────────────────────────────────────────────────────"
    
    if [ $BLAST_EXIT -eq 0 ]; then
        # Move temp file to final location
        mv "$TEMP_OUTPUT" "${OUTPUT_DIR}/${OUTPUT_FILE}"
        
        END_TIME=$(date +%s)
        TOTAL_TIME=$((END_TIME - START_TIME))
        
        echo -e "${GREEN}✓ BLAST search completed successfully${NC}"
        echo -e "   Results: ${BLUE}${OUTPUT_DIR}/${OUTPUT_FILE}${NC}"
        echo "   Total hits: $(wc -l < ${OUTPUT_DIR}/${OUTPUT_FILE})"
        echo "   Total time: $(printf '%02d:%02d:%02d' $((TOTAL_TIME/3600)) $(((TOTAL_TIME%3600)/60)) $((TOTAL_TIME%60)))"
        echo "   Progress log: $PROGRESS_LOG"
        
        echo "Completed: $(date)" >> "$PROGRESS_LOG"
        echo "Total time: ${TOTAL_TIME}s" >> "$PROGRESS_LOG"
    else
        echo -e "${RED}ERROR:${NC} BLAST failed with exit code $BLAST_EXIT"
        exit $BLAST_EXIT
    fi
    
    current_timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "   Completed at: $current_timestamp"
    
    # Clean up temp files
    rm -f "${OUTPUT_DIR}/.blast_results_temp.tsv"
    # remove the peptide_db files
    rm -f "${OUTPUT_DIR}/${DB_NAME}."*
fi


exit 0
