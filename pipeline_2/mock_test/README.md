# Pipeline 2 Testing

Simple testing framework for Pipeline 2 Ks calculation workflow.

## Files

- **`run_test.sh`** - Main test script (synthetic or real data)
- **`make_test_data.py`** - Creates test datasets from clustering results

## Quick Test

```bash
# Fast test with synthetic data (2-3 minutes)
./run_test.sh synthetic

# Test with real data (5-10 minutes)  
./run_test.sh real
```

## What it tests

1. **Step 1**: Extract pairwise combinations from gene families
2. **Step 2**: Align protein sequences with ClustalW2
3. **Step 3**: Back-translate to codon alignments with pal2nal
4. **Step 4**: Calculate Ks/Ka values with PAML yn00
5. **Step 5**: Consolidate and filter results

## Requirements

- Pipeline 1 completed (for real data test)
- ClustalW2, PAML yn00 installed
- Python 3 with pandas

## Output

Test creates `test_data/` directory with:
- Input files (families, sequences)
- Pipeline outputs (alignments, Ks results)
- Log files for each step