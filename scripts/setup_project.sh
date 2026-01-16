#!/bin/bash
# Setup script for Glycine Max comparative genomics project
# Creates necessary directories and downloads initial data
# Date: 2025-11-23

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Glycine Max Comparative Genomics - Project Setup ===${NC}"
echo ""

sudo apt install ncbi-blast+
sudo apt-get install mcl  

# Configuration
RELEASE=41
SPECIES="glycine_max"
BASE_URL="http://ftp.ensemblgenomes.org/pub/release-${RELEASE}/plants"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "Project root: ${PROJECT_ROOT}"
echo "Data source: Ensembl Plants Release ${RELEASE}"
echo ""

# Create directory structure
echo -e "${YELLOW}Creating directory structure...${NC}"

DIRS=(
    "data/genome"
    "data/annotations"
    "data/TE"
    "data/reference"
    "output/statistics"
    "output/ks_values"
    "output/TE_coverage"
    "output/figures"
    "logs"
)

for dir in "${DIRS[@]}"; do
    mkdir -p "${PROJECT_ROOT}/${dir}"
    echo "  ✓ ${dir}"
done

echo ""

# Check if data already exists
GENOME_FILE="${PROJECT_ROOT}/data/genome/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz"
GFF_FILE="${PROJECT_ROOT}/data/annotations/Glycine_max.Glycine_max_v2.1.${RELEASE}.gff3.gz"
PROTEIN_FILE="${PROJECT_ROOT}/data/genome/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz"
CDS_FILE="${PROJECT_ROOT}/data/genome/Glycine_max.Glycine_max_v2.1.cds.all.fa.gz"

# Download genome data (optional - can be large)
read -p "Download genome FASTA file? (~350 MB) [y/N]: " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Downloading genome FASTA...${NC}"
    GENOME_URL="${BASE_URL}/fasta/${SPECIES}/dna/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz"
    wget -O "${GENOME_FILE}" "${GENOME_URL}" || echo -e "${RED}Download failed!${NC}"
fi

# Download protein sequences
read -p "Download protein sequences? (~15 MB) [y/N]: " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Downloading protein sequences...${NC}"
    PROTEIN_URL="${BASE_URL}/fasta/${SPECIES}/pep/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz"
    wget -O "${PROTEIN_FILE}" "${PROTEIN_URL}" || echo -e "${RED}Download failed!${NC}"
fi

# Download CDS sequences
read -p "Download CDS sequences? (~25 MB) [y/N]: " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Downloading CDS sequences...${NC}"
    CDS_URL="${BASE_URL}/fasta/${SPECIES}/cds/Glycine_max.Glycine_max_v2.1.cds.all.fa.gz"
    wget -O "${CDS_FILE}" "${CDS_URL}" || echo -e "${RED}Download failed!${NC}"
fi

# Download GFF annotation (essential)
if [ ! -f "${GFF_FILE}" ]; then
    echo -e "${YELLOW}Downloading GFF3 annotation...${NC}"
    GFF_URL="${BASE_URL}/gff3/${SPECIES}/Glycine_max.Glycine_max_v2.1.${RELEASE}.gff3.gz"
    wget -O "${GFF_FILE}" "${GFF_URL}" || echo -e "${RED}Download failed!${NC}"
else
    echo -e "${GREEN}GFF3 file already exists.${NC}"
fi

# Generate chromosome sizes file
if [ -f "${GENOME_FILE}" ]; then
    echo -e "${YELLOW}Generating chromosome sizes file...${NC}"
    
    # Check if samtools is available
    if command -v samtools &> /dev/null; then
        samtools faidx "${GENOME_FILE}"
        cut -f1,2 "${GENOME_FILE}.fai" > "${PROJECT_ROOT}/data/reference/chromosome_sizes.tsv"
        echo -e "${GREEN}  ✓ Chromosome sizes created${NC}"
    else
        echo -e "${RED}  ✗ samtools not found. Trying alternative method...${NC}"
        
        # Try extracting from GFF
        if [ -f "${GFF_FILE}" ]; then
            zcat "${GFF_FILE}" | grep "^##sequence-region" | \
                awk '{print $2"\t"$4}' > "${PROJECT_ROOT}/data/reference/chromosome_sizes.tsv"
            echo -e "${GREEN}  ✓ Chromosome sizes extracted from GFF${NC}"
        fi
    fi
fi

# Filter main chromosomes
if [ -f "${PROJECT_ROOT}/data/reference/chromosome_sizes.tsv" ]; then
    echo -e "${YELLOW}Filtering main chromosomes...${NC}"
    
    # Glycine max chromosomes: Gm01-Gm20
    grep -E "^Gm[0-9]{2}" "${PROJECT_ROOT}/data/reference/chromosome_sizes.tsv" \
        > "${PROJECT_ROOT}/data/reference/main_chromosomes.tsv" || true
    
    CHR_COUNT=$(wc -l < "${PROJECT_ROOT}/data/reference/main_chromosomes.tsv")
    echo -e "${GREEN}  ✓ Found ${CHR_COUNT} main chromosomes${NC}"
fi

# Extract organellar genes (to exclude)
if [ -f "${GFF_FILE}" ]; then
    echo -e "${YELLOW}Extracting organellar gene list...${NC}"
    
    zcat "${GFF_FILE}" | grep -E "^(Mt|Pt|ChrC|ChrM)" | \
        awk '$3=="gene"' | cut -f9 | \
        grep -o "ID=[^;]*" | sed 's/ID=//' \
        > "${PROJECT_ROOT}/data/reference/organellar_genes.txt" || true
    
    ORG_COUNT=$(wc -l < "${PROJECT_ROOT}/data/reference/organellar_genes.txt" 2>/dev/null || echo 0)
    echo -e "${GREEN}  ✓ Found ${ORG_COUNT} organellar genes to exclude${NC}"
fi

# Create data provenance log
echo -e "${YELLOW}Creating data provenance log...${NC}"

LOG_FILE="${PROJECT_ROOT}/logs/data_provenance.log"
cat > "${LOG_FILE}" <<EOF
# Data Provenance Log
# Glycine Max Comparative Genomics Project

Date: $(date)
Data Source: Ensembl Plants
Release: ${RELEASE}
Species: ${SPECIES}
Genome Version: Glycine_max_v2.1

## Downloaded Files

Base URL: ${BASE_URL}

Genome FASTA:
  ${BASE_URL}/fasta/${SPECIES}/dna/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz

Protein sequences:
  ${BASE_URL}/fasta/${SPECIES}/pep/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz

CDS sequences:
  ${BASE_URL}/fasta/${SPECIES}/cds/Glycine_max.Glycine_max_v2.1.cds.all.fa.gz

GFF3 annotation:
  ${BASE_URL}/gff3/${SPECIES}/Glycine_max.Glycine_max_v2.1.${RELEASE}.gff3.gz

## Reference Files Generated

- chromosome_sizes.tsv: All chromosome/scaffold sizes
- main_chromosomes.tsv: Main chromosomes only (Gm01-Gm20)
- organellar_genes.txt: Mitochondrial and chloroplast genes (to exclude)

## Notes

- All data from same genome release for version consistency
- TE annotation will be downloaded separately (see analysis/TE_analysis/)
- Protein and CDS files needed for pipeline_2 (Ks calculation)

EOF

echo -e "${GREEN}  ✓ Data provenance log created${NC}"

# Summary
echo ""
echo -e "${GREEN}=== Setup Complete ===${NC}"
echo ""
echo "Directory structure created ✓"
echo "Reference files prepared ✓"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "1. Run pipeline_1 if not already done (duplicated gene identification)"
echo "2. Run pipeline_2 if not already done (Ks calculation)"
echo "3. Download TE annotation: cd analysis/TE_analysis && bash 01_download_TE.sh"
echo "4. Start analysis: see analysis/README.md for workflows"
echo ""
echo -e "${YELLOW}Documentation:${NC}"
echo "  - Project plan: docs/project_plan.md"
echo "  - Research questions: docs/research_questions.md"
echo "  - Analysis guide: analysis/README.md"
echo "  - Task list: TODO.md"
echo ""
echo -e "${GREEN}Happy analyzing!${NC}"
