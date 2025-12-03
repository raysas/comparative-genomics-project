# M2 GENIOMHE - Comparative Genomics

## Project Information

**Organism:** *Glycine max* (Soybean)  
**Genome Source:** Ensembl Plants Release 41  
**URL:** http://ftp.ensemblgenomes.org/pub/release-41/plants/fasta/glycine_max/  
**Timeline:** November-December

---

## Research Questions

### Primary Questions
1. What is the extent of gene duplication in the Glycine max genome?
2. What are the ages and mechanisms of these duplications?
3. How are transposable elements distributed across the genome?
4. What is the relationship between TEs and duplicated genes?
5. Do TAGs (tandemly arrayed genes) differ from other duplicates in age and function?

### Specific TAG Questions
- Are TAGs younger than other types of duplications?
- Do younger TAGs have specific functions?
- Is there a link between duplication age and gene function?

### TE-related Questions
- What is the TE density in different genomic regions?
- Are TEs enriched near duplicated genes?
- How do TE distributions correlate with gene density?

---

## Methodology Overview

### Phase 1: Data Acquisition & Preprocessing

#### Tasks
1. **Download genome data** (Ensembl Plants R41)
   - Protein sequences (FASTA)
   - CDS sequences (FASTA)
   - Gene annotations (GFF3)
   - TE annotations (GFF3) - **ensure version match**

2. **Data cleaning**
   - Remove mitochondrial genes
   - Remove chloroplast genes
   - Filter scaffolds (keep only whole chromosomes)
   - Select longest isoform per gene

3. **Create reference files**
   - Gene → Chromosome → Position mapping
   - Gene → Protein → CDS mapping
   - Organellar gene list (for exclusion)

---

### Phase 2: Duplicated Gene Identification (Pipeline 1)

#### Workflow
```
Raw genome → Filter isoforms → BLAST → Coverage → Filter → Network → MCL clustering
```

#### Parameters
- BLAST: E-value ≤ 1e-5
- Identity: ≥ 30%
- Query coverage: ≥ 50%
- Subject coverage: ≥ 50%
- MCL inflation: 1.4

**Outputs:**
- `output/clusters/protein_families_*.tsv` - gene to family mapping
- `output/clusters/protein_families_*.txt` - MCL format families

---

### Phase 3: Ks Calculation (Pipeline 2)

#### Workflow
```
Duplicated pairs → Align proteins → Back-translate → Phylip format → PAML yn00 → Ks values
```

#### Tools
- MUSCLE or MAFFT: protein alignment
- pal2nal: codon alignment
- PAML yn00: Ks/Ka calculation

#### Quality Control
- Remove pairs with Ks > 5 (saturation)
- Remove pairs with alignment length < 100 codons
- Flag pairs with Ka/Ks > 1 (potential non-functional)

**Outputs:**
- Ks values per gene pair (TSV)
- Alignment statistics
- QC report

---

### Phase 4: Statistical Analysis & Data Integration

#### 4.1 Basic Statistics

Compute:
- Total genes in genome
- Total duplicated genes
- Number of families
- Family size distribution
- Genes per chromosome
- Duplicated genes per chromosome

**Outputs:**
- Summary statistics table
- Distribution plots

---

#### 4.2 Duplication Type Classification

**Criteria:**
1. **Tandem (TAG):** Same chromosome, ≤10 genes apart OR ≤100kb distance
2. **Proximal:** Same chromosome, >100kb but <500kb apart
3. **Dispersed:** Different chromosomes OR >500kb apart
4. **WGD:** Ks peak corresponding to known WGD events

**Outputs:**
- Gene annotation file with duplication type
- Classification statistics

---

#### 4.3 Ks Distribution Analysis

**Analysis:**
1. Overall Ks distribution (histogram, density plot)
2. Ks by duplication type
3. Peak identification (mixture models)
4. WGD event detection

**Known WGD events in Glycine max:**
- ~13 MYA (Glycine-specific)
- ~59 MYA (Papilionoideae-wide)

**Conversion:** Age (MYA) = Ks / (2 × λ)  
where λ = synonymous substitution rate (~6.5 × 10⁻⁹ per site per year for legumes)

**Outputs:**
- Ks distribution plots
- Peak analysis results
- Age estimates table

---

### Phase 5: TE Analysis

#### 5.1 TE Annotation Retrieval

- TE annotation GFF
- Version verification

---

#### 5.2 TE Parsing & Classification

**Extract:**
- TE coordinates (chr, start, end, strand)
- TE type/superfamily (LTR, LINE, SINE, DNA, etc.)
- TE length
- TE family

**Outputs:**
- Parsed TE table (TSV)
- TE classification summary

---

#### 5.3 TE Distribution Analysis

**Metrics:**
1. **By type:**
   - Count per superfamily
   - Length distribution per type
   - Genome coverage per type

2. **By chromosome:**
   - TE density per chromosome
   - TE coverage (%) per chromosome

3. **By genomic feature:**
   - Distance to nearest gene
   - TE density in genic vs intergenic regions

**Outputs:**
- TE statistics table
- Distribution plots
- TE density per chromosome plot

---

#### 5.4 TE Coverage Computation (Window-based)

**Algorithm:**
```python
# For each chromosome:
#   1. Define bins/windows (e.g., 100kb)
#   2. For each bin:
#      - Count overlapping TEs
#      - Sum TE coverage (bp)
#      - Calculate density (TE_bp / bin_size)
#   3. Annotate bin features (gene density, duplication density)
```

**Window sizes to test:** 50kb, 100kb, 500kb, 1Mb

**Outputs:**
- TE coverage per bin (TSV)
- Coverage heatmaps
- Correlation with gene density

---

#### 5.5 TE Density Testing

**Comparisons:**
1. TE density near duplicated genes vs singletons
2. TE density near TAGs vs non-TAGs
3. TE density in gene-rich vs gene-poor regions
4. TE density near young vs old duplications

**Statistical tests:**
- Mann-Whitney U test
- Permutation tests for significance

---

### Phase 6: TAG Analysis

#### 6.1 TAG Identification

**Criteria:**
- Same chromosome
- Within 10 genes OR within 100kb
- Both genes in same family

**Output:**
- TAG clusters (TSV)
- TAG vs non-TAG gene list

---

#### 6.2 TAG Orientation Analysis

**Analysis:**
1. Extract gene orientations from GFF
2. For each TAG pair, determine:
   - Same orientation (+/+ or -/-)
   - Opposite orientation (+/-)
3. Test against null hypothesis (random orientation)
4. Statistical test: χ² goodness-of-fit


#### 6.3 TAG Age Analysis

**Comparisons:**
- Ks distribution: TAGs vs WGD vs dispersed
- Mean/median Ks comparison
- Statistical tests (Mann-Whitney, Kruskal-Wallis)

**Hypothesis:** TAGs are younger (lower Ks) than other duplication types

**Outputs:**
- Ks distribution plots by type
- Age summary table

---

#### 6.4 TAG Functional Analysis

**Tools:** clusterProfiler (R) or Panther

**Analysis:**
1. **GO enrichment:**
   - TAGs vs non-TAGs
   - Young TAGs (Ks < 0.5) vs old TAGs
   - Compare enriched GO terms

2. **Age-function correlation:**
   - Bin genes by Ks
   - Test for functional enrichment per bin
   - Identify age-specific functions

**Outputs:**
- GO enrichment tables
- Enrichment plots (dotplots, barplots)
- Age-function analysis results

---

## Technical Specifications

### Software Requirements

#### Core Tools
- **BLAST+** 2.10+: sequence similarity
- **MCL** 14.x: graph clustering
- **MUSCLE** 3.8+ or **MAFFT** 7.x: alignment
- **pal2nal** v14: codon alignment
- **PAML** 4.9: Ks/Ka calculation

---

## Directory Structure
```
analysis/
  ├── duplicated_genes/
  │   ├── TAGs/
  │   └── [statistics scripts]
  ├── Ks/
  │   └── [Ks analysis scripts]
  └── TE_analysis/
      └── [TE scripts and data]

data/
  ├── genome/              # Raw genome files
  ├── annotations/         # GFF files
  ├── TE/                  # TE annotations
  └── reference/           # Gene mappings

output/
  ├── clusters/            # Gene families
  ├── ks_values/           # Ks calculations
  ├── statistics/          # Summary statistics
  ├── TE_coverage/         # TE analysis results
  └── figures/             # All plots
```

---

## References & Resources

### Key Literature on Glycine max
1. Schmutz et al. (2010) "Genome sequence of the palaeopolyploid soybean" *Nature*
2. Known WGD events: ~13 MYA and ~59 MYA
3. Soybean genome duplication and evolution studies

### Data Sources
- Ensembl Plants: http://plants.ensembl.org/
- STRING PPI database: https://string-db.org/
- Panther GO: http://www.pantherdb.org/


---
