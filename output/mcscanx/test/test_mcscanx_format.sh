#!/bin/bash
# Quick test to verify MCScanX format requirements
# Creates minimal test files and tries to run MCScanX

set -euo pipefail

TEST_DIR="test"
mkdir -p "$TEST_DIR"

# Create minimal test GFF (5 genes on 2 chromosomes)
cat > "$TEST_DIR/test.gff" <<'EOF'
gm1	gene1	1000	1500
gm1	gene2	2000	2500
gm1	gene3	3000	3500
gm2	gene4	1000	1500
gm2	gene5	2000	2500
EOF

# Create minimal test BLAST (some hits between genes)
cat > "$TEST_DIR/test.blast" <<'EOF'
gene1	gene2	1e-50	200
gene1	gene4	1e-40	180
gene2	gene3	1e-45	190
gene2	gene5	1e-35	170
gene3	gene5	1e-30	160
EOF

echo "Created test files in $TEST_DIR/"
echo "GFF format:"
head -3 "$TEST_DIR/test.gff"
echo ""
echo "BLAST format:"
head -3 "$TEST_DIR/test.blast"
echo ""

if command -v MCScanX >/dev/null 2>&1; then
    echo "Running MCScanX test..."
    cd "$TEST_DIR"
    MCScanX test 2>&1 | head -20
    
    if [[ -f "test.collinearity" ]]; then
        echo ""
        echo "✓ SUCCESS! MCScanX ran. Collinearity file created:"
        head -10 test.collinearity
    else
        echo ""
        echo "⚠ MCScanX ran but no collinearity found (may need more data)"
    fi
else
    echo "MCScanX not found - skipping test"
fi
