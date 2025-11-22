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
LINK="http://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.0.pep.all.fa.gz"
OUTPUT_DIR='data'
FASTA_FILE="peptides.fa"
FEATURE_OUTPUT_FILE='protein_info.csv'

# -- message on what this script does
cat <<EOF

-- this script downloads peptide data from Ensembl Plants,
   extracts protein information from fasta headers,
   and saves it in a tabular format: 2 outputs
     1) the fasta file with peptides
     2) a tsv file with extracted protein info

EOF

# -- get arguments if provided any
while getopts l:o:i:h flag
do
    case "${flag}" in
        l) LINK=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
        i) FASTA_FILE=${OPTARG};;
        h) echo "Usage: $0 [-l download_link] [-o output_directory] [-i FASTA_FILE]"
           echo "  -l download link for peptide fasta file (default: $LINK)"
           echo "  -o output directory to save  data files in (default: $OUTPUT_DIR)"
           echo "  -i fasta file name to be saved (default: $FASTA_FILE)"
           exit 0;;
    esac
done

# -- define full paths for output files
FASTA_FILE="${OUTPUT_DIR}/${FASTA_FILE}"
FEATURE_OUTPUT_FILE="${OUTPUT_DIR}/${FEATURE_OUTPUT_FILE}"

# -- set up logging
LOG_DIR="logs/pipeline"
if [ ! -d "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi
LOG_FILE="${LOG_DIR}/$(basename "$0" .sh).log"
exec > >(tee -i "$LOG_FILE") 2>&1
echo "Command: $0 $*"

# -- start message
echo "-- downloading peptide data from Ensembl Plants:"
echo "   LINK  : $LINK"
echo "   OUTPUT PATH: $FASTA_FILE"
echo "   FEATURE OUTPUT PATH: $FEATURE_OUTPUT_FILE"

# -- ensure output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
    echo "-- created output directory: $OUTPUT_DIR"
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

# -- 1) if fasta already exists, skip this step
if [ ! -f "$FASTA_FILE" ]; then
    echo "-- fasta file does not exist, proceeding to download."
    date_now=$(date +"%Y-%m-%d %T")
    echo "-- [$date_now] downloading data from $LINK"

    curl $LINK -o "${FASTA_FILE}.gz"
    gunzip "${FASTA_FILE}.gz"

    # -- return the file path
    echo "-- data downloaded to $FASTA_FILE"
else
    echo "-- fasta file already exists at $FASTA_FILE, skipping download."
fi

# -- 2) extract protein information from fasta headers
#    from these headers need to extract all info >KRH46721 pep chromosome:Glycine_max_v2.0:8:46595288:46603445:-1 gene:GLYMA_08G352800 transcript:KRH46721 gene_biotype:protein_coding transcript_biotype:protein_coding description:hypothetical protein

echo " -- now extracting protein information from fasta headers:"
echo "peptide_id,gene_id,transcript_id,genome,chromosome,start_pos,end_pos,strand,description,length" > "$FEATURE_OUTPUT_FILE"

# -- extract to metadata file
grep ">" "$FASTA_FILE" | sed 's/>//g' | awk '{ printf "%s %s %s %s ", $1, $3, $4, $5; for(i=8;i<=NF;i++) {printf "%s%s", $i, (i<NF?OFS:ORS) }}' | awk '{
    match($0, /chromosome:([^ ]+)/, chr);
    split(chr[1], a, ":"); 

    match($0, /gene:([^ ]+)/, gene);
    match($0, /transcript:([^ ]+)/, transcript);
    match($0, /description:(.*)$/, desc);  # Capture EVERYTHING after description:

    print $1 "," gene[1] "," transcript[1] "," a[1] "," a[2] "," a[3] "," a[4] "," a[5] "," desc[1];
}' >> "${OUTPUT_DIR}/protein_metadata.csv"

# -- extract sequence lengths and append to metadata file
awk '/^>/ {if(seq){print length(seq)}; seq=""; next} {seq=seq $0} END{print length(seq)}' "$FASTA_FILE" \
> "${OUTPUT_DIR}/seq_lengths.csv"

paste -d, "${OUTPUT_DIR}/protein_metadata.csv" "${OUTPUT_DIR}/seq_lengths.csv" >> "$FEATURE_OUTPUT_FILE"

# -------------------------------------------------------------------------
# -- end message
echo "-- protein info extracted to $FEATURE_OUTPUT_FILE"
echo "-- can find them in this form: "
head -n 1 "$FEATURE_OUTPUT_FILE"
exit 0

