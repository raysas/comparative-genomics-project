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
cat <<EOF
-- this script:
  1) creates a BLAST database from a given peptide fasta file
  2) runs a BLASTP search of the peptides against themselves
EOF

# -- default parameters
INPUT_FILE=''
OUTPUT_DIR=''
AUTO_DETECT=true

DB_NAME='peptide_db'
OUTPUT_FILE='blast_results.tsv'

# -- make a flag to not keep blast_output directory across runs and make a new one each time
force_new=false

# -- get arguments
while getopts "i:o:fh" flag; do
    case "${flag}" in
        i) INPUT_FILE="${OPTARG}" ;;
        o) OUTPUT_DIR="${OPTARG}" ;;
        f) force_new=true ;;
        h)
            cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -i FILE       Input FASTA file (peptides_longest.fa)
  -o DIR        Output directory for BLAST results
  -f            Force new output directory (rename existing)
  -h            Show this help

AUTO-DETECTION:
  If no options specified, automatically detects species directory:
    Input:  data/{species}/peptides_longest.fa
    Output: output/{species}/blast_output

EXAMPLES:
  # Auto-detect (recommended)
  $0

  # Explicit paths
  $0 -i data/glycine_max/peptides_longest.fa -o output/glycine_max/blast_output

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
            
            [ -z "$INPUT_FILE" ] && INPUT_FILE="${SPECIES_DIR}/peptides_longest.fa"
            [ -z "$OUTPUT_DIR" ] && OUTPUT_DIR="output/${SPECIES_NAME}/blast_output"
            
            echo "   Detected species: $SPECIES_NAME"
        else
            echo "ERROR: Cannot auto-detect. Provide both -i and -o"
            exit 1
        fi
    fi
fi

# Auto-detect species directory if not specified
if [ "$AUTO_DETECT" = true ]; then
    echo "-- Auto-detecting species directory..."
    
    SPECIES_DIRS=(data/*/)
    
    if [ ${#SPECIES_DIRS[@]} -eq 1 ] && [ -d "${SPECIES_DIRS[0]}" ]; then
        SPECIES_DIR="${SPECIES_DIRS[0]%/}"
        SPECIES_NAME=$(basename "$SPECIES_DIR")
        
        INPUT_FILE="${SPECIES_DIR}/peptides_longest.fa"
        OUTPUT_DIR="output/${SPECIES_NAME}/blast_output"
        
        echo "   Detected species: $SPECIES_NAME"
        echo "   Input:  $INPUT_FILE"
        echo "   Output: $OUTPUT_DIR"
    elif [ ${#SPECIES_DIRS[@]} -gt 1 ]; then
        echo "ERROR: Multiple species directories found in data/"
        echo "       Specify files explicitly with -i and -o"
        echo "       Found: ${SPECIES_DIRS[*]}"
        exit 1
    else
        echo "ERROR: No species directories found in data/"
        echo "       Run 0_extract_data.sh and 1_filter_isoforms.sh first"
        exit 1
    fi
fi


# -- name log file based on script name
LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

echo "-- Parameters:"
echo "   INPUT FILE : $INPUT_FILE"
echo "   OUTPUT DIR : $OUTPUT_DIR"
echo "   FORCE NEW  : $force_new"

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
    echo "-- preparing to create new output directory: $OUTPUT_DIR"
    leaf_dir=$(basename "$OUTPUT_DIR")
    parent_dir=$(dirname "$OUTPUT_DIR")
    if [ -d "$OUTPUT_DIR" ]; then
        # -- move leaf to leaf_old_numbered
        n=1
        while [ -d "${parent_dir}/${leaf_dir}_old_${n}" ]; do
            ((n++))
        done
        mv "$OUTPUT_DIR" "${parent_dir}/${leaf_dir}_old_${n}"
        echo "-- ! moved existing directory $OUTPUT_DIR to ${parent_dir}/${leaf_dir}_old_${n}"
    fi

    # -- create output directory
    mkdir -p "$OUTPUT_DIR"
    echo "-- created output directory: $OUTPUT_DIR"

    # -- make blast db
    if ! command -v makeblastdb &> /dev/null
    then
        echo "-- error: makeblastdb could not be found: install BLAST+ tools"
        exit 2
    fi

    makeblastdb -in "$INPUT_FILE" -dbtype prot -out "${OUTPUT_DIR}/${DB_NAME}"
    echo "-- created BLAST database from $INPUT_FILE"
    
    # Count total sequences for progress estimation
    TOTAL_SEQS=$(grep -c "^>" "$INPUT_FILE")
    echo "-- Total sequences to process: $TOTAL_SEQS"

    # -- run blastp with optimizations for speed
    if ! command -v blastp &> /dev/null
    then
        echo "-- error: blastp could not be found: install BLAST+ tools"
        exit 3
    fi
    
    # Detect available CPU cores
    NCORES=$(nproc 2>/dev/null || echo 4)
    # Use 75% of cores for BLAST
    BLAST_THREADS=$((NCORES * 3 / 4))
    [ $BLAST_THREADS -lt 1 ] && BLAST_THREADS=1
    
    echo "-- Running BLASTP with $BLAST_THREADS threads (detected $NCORES cores)"
    echo "   This may take a while for large genomes..."
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
           -evalue 1e-5 \
           -max_target_seqs 10000 \
           -num_threads $BLAST_THREADS \
           -word_size 6 \
           -threshold 21 \
           -comp_based_stats 0 \
           -seg no &
    
    BLAST_PID=$!
    
    echo "BLAST process started (PID: $BLAST_PID)"
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
        
        echo "-- BLAST search completed successfully"
        echo "   Results: ${OUTPUT_DIR}/${OUTPUT_FILE}"
        echo "   Total hits: $(wc -l < ${OUTPUT_DIR}/${OUTPUT_FILE})"
        echo "   Total time: $(printf '%02d:%02d:%02d' $((TOTAL_TIME/3600)) $(((TOTAL_TIME%3600)/60)) $((TOTAL_TIME%60)))"
        echo "   Progress log: $PROGRESS_LOG"
        
        echo "Completed: $(date)" >> "$PROGRESS_LOG"
        echo "Total time: ${TOTAL_TIME}s" >> "$PROGRESS_LOG"
    else
        echo "ERROR: BLAST failed with exit code $BLAST_EXIT"
        exit $BLAST_EXIT
    fi
    
    current_timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "   Completed at: $current_timestamp"
    
    # Clean up temp files
    rm -f "${OUTPUT_DIR}/.blast_results_temp.tsv"
fi


exit 0
