#!/bin/bash

# Minimal Pipeline 2 Test Script
# Tests all 5 pipeline steps with real or synthetic data

set -e

echo "=== Pipeline 2 Test ==="

# Test mode: synthetic (fast) or real (uses actual clustering data)
MODE=${1:-synthetic}

if [ "$MODE" = "real" ]; then
    echo "🧬 Testing with REAL data"
    
    # Check prerequisites
    FULL_FILE="../../output/clusters/protein_families_filtered_blast_results_id30_qcov50_scov50_wcol12_network.tsv"
    if [ ! -f "$FULL_FILE" ]; then
        echo "❌ Missing clustering results. Run pipeline_1 first."
        exit 1
    fi
    
    # Create test dataset
    echo "📊 Creating test dataset..."
    python3 make_test_data.py "$FULL_FILE" "test_data/families.tsv" "small"
    
    # Extract protein sequences for test genes (optimized)
    echo "🧬 Extracting sequences..."
    
    # Create gene ID lookup for faster extraction
    tail -n +2 test_data/families.tsv | cut -f1 | sort > test_data/gene_ids.txt
    
    # Extract proteins (single pass)
    awk 'BEGIN {
        while((getline id < "test_data/gene_ids.txt") > 0) genes[id]=1; close("test_data/gene_ids.txt")
    }
    /^>/ {
        gene_id = $1; gsub(/^>/, "", gene_id)
        if(gene_id in genes) flag=1; else flag=0
    }
    flag' "../../data/glycine_max/peptides_longest.fa" > test_data/proteins.fa
    
    # Extract CDS if available (single pass)
    if [ -f "../../data/cds.fa" ]; then
        awk 'BEGIN {
            while((getline id < "test_data/gene_ids.txt") > 0) genes[id]=1; close("test_data/gene_ids.txt")
        }
        /^>/ {
            gene_id = $1; gsub(/^>/, "", gene_id)
            if(gene_id in genes) flag=1; else flag=0
        }
        flag' "../../data/cds.fa" > test_data/cds.fa
    else
        echo "⚠ No CDS file - will test only protein alignment steps"
        touch test_data/cds.fa
    fi
    
    # Clean up temporary file
    rm -f test_data/gene_ids.txt
    
    # Set paths for pipeline (using absolute paths)
    export FAMILIES_FILE="$(pwd)/test_data/families.tsv"
    export PROTEIN_FILE="$(pwd)/test_data/proteins.fa" 
    export CDS_FILE="$(pwd)/test_data/cds.fa"
    export OUTPUT_DIR="$(pwd)/test_data/output"
    export SPECIES="glycine_max"
    
else
    echo "⚡ Testing with SYNTHETIC data (fast)"
    
    # Create minimal synthetic data
    mkdir -p test_data
    
    # Synthetic families (2 families, 2-3 genes each)
    cat > test_data/families.tsv << 'EOF'
geneName	family
TEST001	1
TEST002	1
TEST003	2
TEST004	2
TEST005	2
EOF
    
    # Synthetic proteins (similar sequences for alignment)
    cat > test_data/proteins.fa << 'EOF'
>TEST001
MKAILVVLLYTFATAAQAEDRPSTRDMACPLIIQVLSLLARAKAKKKKKEDNGKGKGQQLHKPVGMQHFRIVDGAATIYSVAQVTKDDMDDDKTKDDPK
>TEST002
MKAILVVLLYTFATAAQAEDRPSTRDMACPLVIQVLSLLARAKAKKKKKEDNGKGKGQQLHKPVGMQHFRIVDGAATIYSVAQVTKDDMDDDKTKDDPK
>TEST003
MHNLGPIVETLGIAGNRNFRQALKRFSQSGKASILQTVFSDGEGKQDWWKAAGSQTVIEKATGKDYDFVPITAMQWIAGHKLRAWDIDEFFYESNFK
>TEST004
MHNLGPIVETLGIAGNRNFRQALKRFSQSGKASILQTVFSDGEGKQDWWKAAGSQTVIEKATGKDYDFVPITAMQWIAGHKLRAWDIDEFFYESNFK
>TEST005
MHNLGPIVETLGIAGNRNFRQALKRFSQSGKASILQTVFSDGEGKQDWWKAAGSQTVIEKATGKDYDFVPITAMQWIAGHKLRAWDIDEFFYESNFQ
EOF
    
    # Synthetic CDS (matching proteins)
    cat > test_data/cds.fa << 'EOF'
>TEST001
ATGAAAGCTATTCTGGTAGTGTTGTTGTACACTTTCGCTACAGCTGCTCAAGCTGAAGACCGTCCTTCAACTCGTGACATGGCCTGCCCTTTAATTATCCAAGTGCTTTCGTTGCTGGCTCGTGCAAAAGCAAAGAAAAAGAAAAAGAAAGACAACGGCAAAGGCAAAGGCCAACAATTACATAAACCAGTAGGCATGCAACACTTTAGAATAGTGGACGGCGCAGCTACAATCTATAGCGTAGCACAAGTGACTAAAGACGATATGGACGACGACAAAACTAAAGACGATCCGAAA
>TEST002
ATGAAAGCTATTCTGGTAGTGTTGTTGTACACTTTCGCTACAGCTGCTCAAGCTGAAGACCGTCCTTCAACTCGTGACATGGCCTGCCCTTTGGTTATCCAAGTGCTTTCGTTGCTGGCTCGTGCAAAAGCAAAGAAAAAGAAAAAGAAAGACAACGGCAAAGGCAAAGGCCAACAATTACATAAACCAGTAGGCATGCAACACTTTAGAATAGTGGACGGCGCAGCTACAATCTATAGCGTAGCACAAGTGACTAAAGACGATATGGACGACGACAAAACTAAAGACGATCCGAAA
>TEST003
ATGCACAACCTGGGCCCAATCGTGGAAACGCTGGGCATTGCGGGCAACCGCAACTTCCGCCAGGCGCTGAAACGCTTCAGCCAGTCCGGCAAAGCGTCCATCCTGCAAACCGTGTTTAGCGACGGCGAAGGCAAACAGGACTGGTGGAAAGCGGCGGGCAGCCAAACCGTGATCGAGAAAGCGACCGGCAAAGACTACGACTTCGTGCCAATCACCGCGATGCAGTGGATTGCGGGCCACAAACTGCGCGCGTGGGATATTGACGAGTTCTTCTACGAGAGCAACTTCAAA
>TEST004
ATGCACAACCTGGGCCCAATCGTGGAAACGCTGGGCATTGCGGGCAACCGCAACTTCCGCCAGGCGCTGAAACGCTTCAGCCAGTCCGGCAAAGCGTCCATCCTGCAAACCGTGTTTAGCGACGGCGAAGGCAAACAGGACTGGTGGAAAGCGGCGGGCAGCCAAACCGTGATCGAGAAAGCGACCGGCAAAGACTACGACTTCGTGCCAATCACCGCGATGCAGTGGATTGCGGGCCACAAACTGCGCGCGTGGGATATTGACGAGTTCTTCTACGAGAGCAACTTCAAA
>TEST005
ATGCACAACCTGGGCCCAATCGTGGAAACGCTGGGCATTGCGGGCAACCGCAACTTCCGCCAGGCGCTGAAACGCTTCAGCCAGTCCGGCAAAGCGTCCATCCTGCAAACCGTGTTTAGCGACGGCGAAGGCAAACAGGACTGGTGGAAAGCGGCGGGCAGCCAAACCGTGATCGAGAAAGCGACCGGCAAAGACTACGACTTCGTGCCAATCACCGCGATGCAGTGGATTGCGGGCCACAAACTGCGCGCGTGGGATATTGACGAGTTCTTCTACGAGAGCAACTTCCAA
EOF
    
    # Set paths for pipeline (using absolute paths)
    export FAMILIES_FILE="$(pwd)/test_data/families.tsv"
    export PROTEIN_FILE="$(pwd)/test_data/proteins.fa"
    export CDS_FILE="$(pwd)/test_data/cds.fa"
    export OUTPUT_DIR="$(pwd)/test_data/output"
    export SPECIES="test_species"
fi

# Test all pipeline steps with explicit arguments
echo ""
echo "🚀 Running pipeline steps..."

echo "▶ Step 1: Prepare pairs"
if timeout 300 bash ../1_prepare_pairs.sh -i "$FAMILIES_FILE" -p "$PROTEIN_FILE" -o "$OUTPUT_DIR/$SPECIES/pairs" >/dev/null 2>&1; then
    echo "  ✅ Completed"
else
    echo "  ❌ Failed - check ../logs/pipeline/1_prepare_pairs.log"
    if [ -f "../logs/pipeline/1_prepare_pairs.log" ]; then
        tail -5 "../logs/pipeline/1_prepare_pairs.log" | sed 's/^/     /'
    fi
    exit 1
fi

echo "▶ Step 2: Align proteins"  
if timeout 300 bash ../2_align_proteins.sh -i "$OUTPUT_DIR/$SPECIES/pairs" -o "$OUTPUT_DIR/$SPECIES/alignments" >/dev/null 2>&1; then
    echo "  ✅ Completed"
else
    echo "  ❌ Failed - check ../logs/pipeline/2_align_proteins.log"
    if [ -f "../logs/pipeline/2_align_proteins.log" ]; then
        tail -5 "../logs/pipeline/2_align_proteins.log" | sed 's/^/     /'
    fi
    exit 1
fi

echo "▶ Step 3: Back-translate"
if timeout 300 bash ../3_backtranslate.sh -i "$OUTPUT_DIR/$SPECIES/alignments" -o "$OUTPUT_DIR/$SPECIES/codon_alignments" -c "$CDS_FILE" >/dev/null 2>&1; then
    echo "  ✅ Completed"
else
    echo "  ❌ Failed - check ../logs/pipeline/3_backtranslate.log"
    if [ -f "../logs/pipeline/3_backtranslate.log" ]; then
        tail -5 "../logs/pipeline/3_backtranslate.log" | sed 's/^/     /'
    fi
    exit 1
fi

echo "▶ Step 4: Calculate Ks"
if timeout 300 bash ../4_calculate_ks.sh -i "$OUTPUT_DIR/$SPECIES/codon_alignments" -o "$OUTPUT_DIR/$SPECIES/ks_results" >/dev/null 2>&1; then
    echo "  ✅ Completed"
else
    echo "  ❌ Failed - check ../logs/pipeline/4_calculate_ks.log"
    if [ -f "../logs/pipeline/4_calculate_ks.log" ]; then
        tail -5 "../logs/pipeline/4_calculate_ks.log" | sed 's/^/     /'
    fi
    exit 1
fi

echo "▶ Step 5: Consolidate results"
if timeout 300 bash ../5_consolidate_results.sh -i "$OUTPUT_DIR/$SPECIES/ks_results" -o "$OUTPUT_DIR/$SPECIES/ks_summary.tsv" >/dev/null 2>&1; then
    echo "  ✅ Completed"
else
    echo "  ❌ Failed - check ../logs/pipeline/5_consolidate_results.log"
    if [ -f "../logs/pipeline/5_consolidate_results.log" ]; then
        tail -5 "../logs/pipeline/5_consolidate_results.log" | sed 's/^/     /'
    fi
    exit 1
fi

# Show results
echo ""
echo "🎉 All steps completed successfully!"

if [ -f "test_data/output/$SPECIES/ks_summary.tsv" ]; then
    result_count=$(tail -n +2 "test_data/output/$SPECIES/ks_summary.tsv" | wc -l)
    echo "📊 Generated $result_count Ks results"
    echo ""
    echo "Sample results:"
    head -3 "test_data/output/$SPECIES/ks_summary.tsv" | sed 's/^/  /'
fi

echo ""
echo "📁 Output files:"
find test_data/output -name "*.tsv" -o -name "*.aln" -o -name "*.out" | head -5 | sed 's/^/  /'

# Cleanup option
echo ""
read -p "Clean up test data? (Y/n): " -r
if [[ ! $REPLY =~ ^[Nn]$ ]]; then
    rm -rf test_data/
    echo "✅ Cleaned up"
else
    echo "📂 Test data preserved in test_data/"
fi