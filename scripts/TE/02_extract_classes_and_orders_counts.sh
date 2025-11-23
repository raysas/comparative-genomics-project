#!/bin/bash
# ------------------------------------------------------------------
# Script Name: 02_extract_classes_and_orders_counts.sh
# Description: Extracts TE classification from APTE GFF3.
#              1. Fixes Class I/II overlap logic.
#              2. Counts occurrences.
#              3. SORTS by abundance (Highest count at top).
#
# Input:       ../data/TEAnnotationFinal.gff3
# Output:      ../output/TE_classes_and_orders_counts.txt
# ------------------------------------------------------------------

# 1. Define Paths
INPUT_GFF="../data/TEAnnotationFinal.gff3"
OUTPUT_DIR="../output"
OUTPUT_FILE="${OUTPUT_DIR}/TE_classes_and_orders_counts.txt"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

echo "Starting extraction from $INPUT_GFF..."

# 2. Processing Pipeline
grep -v "#" "$INPUT_GFF" | cut -f3 | \
awk -F'/' '{
    # Get the first part (e.g., "Class II subclass 1")
    class_string = $1

    # --- CLASSIFICATION LOGIC ---
    # CRITICAL: Check "Class II" FIRST. 
    # Because regex "/Class I/" matches inside "Class II" (Class I...I).
    if (class_string ~ /Class II/) { 
        group="Class_II" 
    }
    else if (class_string ~ /Class I/) { 
        group="Class_I" 
    }
    else { 
        group="Unknown" 
    }

    # --- ORDER LOGIC ---
    order = $2
    
    # Handle clean up
    if (order == "" || order == "-" || order == " ") { 
        order="Unknown" 
    }

    # Print temp format: Group <tab> Order
    print group "\t" order
}' | \
sort | uniq -c | \
sort -nr | \
awk '{print $1, $2, $3}' > "$OUTPUT_FILE"

# ^^^ EXPLANATION OF THE END OF THE PIPELINE ^^^
# 1. sort       -> Groups identical lines so uniq can count them.
# 2. uniq -c    -> Counts them (Output: "   50 Class_I LTR").
# 3. sort -nr   -> SORTS NUMERICALLY REVERSE (Biggest number first).
# 4. awk        -> Cleans up formatting (removes leading spaces).

# 3. Validation
echo "------------------------------------------------"
echo "Extraction Complete. Data is sorted by abundance."
echo "TE Orders found:"
cat "$OUTPUT_FILE"
echo "------------------------------------------------"