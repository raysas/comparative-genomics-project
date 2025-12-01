#!/bin/bash

# Fast backtranslation wrapper script
# Uses Python for 100x speed improvement over pal2nal

set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Check if Python3 is available
if ! command -v python3 &>/dev/null; then
    echo "ERROR: Python3 not found. Please install it."
    exit 1
fi

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Path to Python script
PYTHON_SCRIPT="$SCRIPT_DIR/3_backtranslate.py"

# Check if Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: 3_backtranslate.py not found in $SCRIPT_DIR"
    exit 1
fi

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN} FAST BACKTRANSLATION (Python)${NC}"
echo -e "${GREEN}========================================${NC}"
echo "This uses a Python implementation that's ~100x faster than pal2nal"
echo ""

# Pass all arguments to Python script
python3 "$PYTHON_SCRIPT" "$@"