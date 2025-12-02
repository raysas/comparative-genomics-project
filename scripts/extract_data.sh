#!/bin/bash

# --------------------------------------------------------------------
# -- What this script does:
#    1) downloads peptide data from Ensembl Plants
#    2) extracts protein information from fasta headers
#       and saves it in a tabular format: 2 outputs
#         a- the fasta file with peptides
#         b- a tsv file with extracted protein info
# -- temporary files created during processing will be stored in the output directory
#    1. metadata.csv
#    2. seq_lengths.csv
# -- Usage:
#    bash ./pipeline_1/0_extract_data.sh [-l download_link] [-o output_directory] [-i FASTA_FILE]
# -- default (without params) equivalent to:
#    bash ./pipeline_1/0_extract_data.sh -l "http://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.0.pep.all.fa.gz" -o "data" -i "peptides.fa"
# --------------------------------------------------------------------

###########################################################################

# ---------------------------------------------------------------------
# prepare variables, get arguments if any, set up logging
# ---------------------------------------------------------------------

# -- default parameters
# Ensembl configuration
SPECIES="glycine_max"
RELEASE="41"
GENOME_VERSION="Glycine_max_v2.0"  # Full genome assembly name as it appears in filename
DATA_TYPE="pep"  # pep for peptides, cds for coding sequences

# Output configuration
OUTPUT_DIR='data'
FASTA_FILE="peptides.fa"
FEATURE_OUTPUT_FILE='protein_info.csv'
SKIP_DOWNLOAD=false

# -- message on what this script does
cat <<EOF

-- this script downloads peptide data from Ensembl Plants,
   extracts protein information from fasta headers,
   and saves it in a tabular format: 2 outputs
     1) the fasta file with peptides
     2) a tsv file with extracted protein info

EOF

# -- get arguments if provided any
while getopts S:r:v:t:o:i:sh flag
do
    case "${flag}" in
        S) SPECIES=${OPTARG};;
        r) RELEASE=${OPTARG};;
        v) GENOME_VERSION=${OPTARG};;
        t) DATA_TYPE=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        i) FASTA_FILE=${OPTARG};;
        s) SKIP_DOWNLOAD=true;;
        h) cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  -S SPECIES        Species name (default: $SPECIES)
  -r RELEASE        Ensembl release number (default: $RELEASE)
  -v VERSION        Genome version (default: $GENOME_VERSION)
  -t TYPE           Data type: pep or cds (default: $DATA_TYPE)
  -o DIR            Output directory (default: $OUTPUT_DIR)
  -i FILE           Output fasta filename (default: $FASTA_FILE)
  -s                Skip download (use existing file)
  -h                Show this help

Example:
    $0 -S oryza_sativa -r 41 -v 1 -t pep

URL FORMAT:
  http://ftp.ensemblgenomes.org/pub/release-{RELEASE}/plants/fasta/{SPECIES}/{TYPE}/{Species}.{GENOME_VERSION}.{TYPE}.all.fa.gz
  Note: Species name is capitalized in filename (e.g., Glycine_max.Glycine_max_v2.0.pep.all.fa.gz)

EOF
           exit 0;;
    esac
done

echo "-- Constructed download URL: $LINK"

# -- List files in remote FTP directory and find the one ending with pep.all.fa.gz
FTP_DIR="http://ftp.ensemblgenomes.org/pub/release-${RELEASE}/plants/fasta/${SPECIES}/${DATA_TYPE}/"
echo "-- Querying remote FTP directory: $FTP_DIR"
PEP_FILE=$(curl -s "$FTP_DIR" | grep -oP 'href="\K[^"]+pep\.all\.fa\.gz' | head -1)
if [ -z "$PEP_FILE" ]; then
    echo "ERROR: Could not find a file ending with pep.all.fa.gz in $FTP_DIR"
    exit 1
fi
LINK="${FTP_DIR}${PEP_FILE}"
echo "-- Found remote peptide file: $PEP_FILE"
echo "-- Download URL: $LINK"

# -- Create species-specific subdirectory within OUTPUT_DIR
OUTPUT_DIR="${OUTPUT_DIR}/${SPECIES}"

# -- define full paths for output files
FASTA_FILE="${OUTPUT_DIR}/${FASTA_FILE}"
FEATURE_OUTPUT_FILE="${OUTPUT_DIR}/processed/${FEATURE_OUTPUT_FILE}"

# -- set up logging
# LOG_DIR="logs/pipeline"
# if [ ! -d "$LOG_DIR" ]; then
#     mkdir -p "$LOG_DIR"
# fi
# LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
# exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

# -- start message
echo "-- Downloading peptide data from Ensembl Plants:"
echo "   Species:        $SPECIES"
echo "   Release:        $RELEASE"
echo "   Genome version: $GENOME_VERSION"
echo "   Data type:      $DATA_TYPE"
echo "   Download URL:   $LINK"
echo "   Output path:    $FASTA_FILE"
echo "   Features path:  $FEATURE_OUTPUT_FILE"
echo ""

# -- ensure species-specific output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "${OUTPUT_DIR}/processed"
    echo "-- created species directory: $OUTPUT_DIR"
fi

# -------------------------------------------------------------------------
#############################################################################
# -------------------------------------------------------------------------

# -- steps to follow:
# 1) download fasta file if not already present
# 2) extract protein information from fasta and save to tsv file
#   a) features from fasta headers => metadata.csv
#   b) sequence lengths => seq_lengths.csv
# -------------------------------------------------------------------------

# -- 1) Download or use existing fasta file
if [ "$SKIP_DOWNLOAD" = true ]; then
    echo "-- skip download flag set, using existing fasta file"
    if [ ! -f "$FASTA_FILE" ]; then
        echo "ERROR: Skip download requested but $FASTA_FILE does not exist"
        exit 1
    fi
elif [ -f "$FASTA_FILE" ]; then
    echo "-- fasta file already exists at $FASTA_FILE, skipping download."
    echo "   Use -s flag explicitly or remove file to re-download"
else
    echo "-- fasta file does not exist, proceeding to download."
    date_now=$(date +"%Y-%m-%d %T")
    echo "-- [$date_now] downloading data from $LINK"

    # Download with progress bar and resume capability
    curl -L --progress-bar --continue-at - "$LINK" -o "${FASTA_FILE}.gz"
    
    # Check download success
    if [ $? -ne 0 ]; then
        echo "ERROR: Download failed"
        exit 1
    fi
    
    echo "-- decompressing file..."
    gunzip -f "${FASTA_FILE}.gz"
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Decompression failed"
        exit 1
    fi

    # -- return the file path
    echo "-- data downloaded to $FASTA_FILE"
fi

# -- 2) extract protein information from fasta headers
#    from these headers need to extract all info >KRH46721 pep chromosome:Glycine_max_v2.0:8:46595288:46603445:-1 gene:GLYMA_08G352800 transcript:KRH46721 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein

# Skip if feature file already exists and is newer than fasta
if [ -f "$FEATURE_OUTPUT_FILE" ] && [ "$FEATURE_OUTPUT_FILE" -nt "$FASTA_FILE" ]; then
    echo "-- protein info file already up-to-date: $FEATURE_OUTPUT_FILE"
    echo "   Delete it to regenerate"
else
    echo " -- extracting protein information from fasta headers..."
    echo "peptide_id,gene_id,transcript_id,genome,chromosome,start_pos,end_pos,strand,description,length" > "$FEATURE_OUTPUT_FILE"

    # -- extract metadata and sequence lengths in one pass for efficiency
    awk '
    BEGIN {
        # Print header to metadata file
        metadata_file = "'"${OUTPUT_DIR}/processed/protein_metadata.csv"'"
    }
    /^>/ {
        # Save previous sequence length if exists
        if (seq_length > 0) {
            print seq_length
        }
        seq_length = 0
        
        # Parse header
        header = substr($0, 2)  # Remove >
        
        # Extract fields using match
        match(header, /^([^ ]+)/, peptide)
        match(header, /chromosome:([^ ]+)/, chr)
        match(header, /gene:([^ ]+)/, gene)
        match(header, /transcript:([^ ]+)/, transcript)
        match(header, /description:(.*)$/, desc)
        
        # Split chromosome info
        split(chr[1], chr_parts, ":")
        
        # Print to metadata
        print peptide[1] "," gene[1] "," transcript[1] "," chr_parts[1] "," chr_parts[2] "," chr_parts[3] "," chr_parts[4] "," chr_parts[5] "," desc[1] > metadata_file
        next
    }
    {
        # Accumulate sequence length
        seq_length += length($0)
    }
    END {
        # Print final sequence length
        if (seq_length > 0) {
            print seq_length
        }
    }
    ' "$FASTA_FILE" > "${OUTPUT_DIR}/processed/seq_lengths.csv"

    # -- combine metadata and lengths
    paste -d, "${OUTPUT_DIR}/processed/protein_metadata.csv" "${OUTPUT_DIR}/processed/seq_lengths.csv" >> "$FEATURE_OUTPUT_FILE"
    
    echo "-- protein info extracted to $FEATURE_OUTPUT_FILE"
fi

# -------------------------------------------------------------------------
# -- end message and summary
echo ""
echo "========================================="
echo "EXTRACTION COMPLETE"
echo "========================================="
echo "FASTA file:    $FASTA_FILE"
echo "Protein info:  $FEATURE_OUTPUT_FILE"
echo ""
echo "Format:"
head -n 1 "$FEATURE_OUTPUT_FILE"
echo ""
echo "Total proteins: $(tail -n +2 "$FEATURE_OUTPUT_FILE" | wc -l)"
echo "========================================="
exit 0

