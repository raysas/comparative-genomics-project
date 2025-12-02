#!/usr/bin/env python
import pandas as pd
import sys

input_file = "output/blast_filtered/filtered_blast_results_id50_qcov70_scov70.tsv"
output_file = "test/per_protein_top5.tsv"

try:
	df = pd.read_csv(input_file, sep=r"\s+", engine="python")
except Exception:
	df = pd.read_csv(input_file, sep="\t")

# normalize column names (strip whitespace)
df.columns = [c.strip() for c in df.columns]

if 'bitscore' not in df.columns:
	candidates = [c for c in df.columns if 'bitscore' in c.lower()]
	if candidates:
		bits_col = candidates[0]
		df = df.rename(columns={bits_col: 'bitscore'})
	else:
		raise KeyError("bitscore column not found in input file columns: %r" % (df.columns.tolist(),))

df_sorted = df.sort_values(by=["qseqid", "bitscore"], ascending=[True, False])
df_top5 = df_sorted.groupby("qseqid").head(5)
df_top5.to_csv(output_file, sep="\t", index=False)

print(f"-- saved top 5 hits per protein to: {output_file}")
