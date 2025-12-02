#!./usr/bin/env python3

# -- to be ran from main project directory
# -- generates diagnostic plots from BLAST results
# -- results in figures/blast/

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================
# -- parameters
# ==========================
blast_file = "output/blast_output/blast_results_with_coverage.tsv"  # your BLAST file
output_prefix = "blast_viz"

# ==========================
# -- load BLAST results
# ==========================
df = pd.read_csv(blast_file, sep='\s+')
df.columns = df.columns.str.strip()  # remove spaces from column names

print('Columns in BLAST file:', df.columns.tolist() )

# Compute minimum coverage between query and subject
df["min_cov"] = df[["qcov", "scov"]].min(axis=1)

print(f"Total hits: {len(df)}")
print(f"Mean percent identity: {df['pident'].mean():.2f}")
print(f"Mean min coverage: {df['min_cov'].mean():.2f}")

# ==========================
# -- plot distributions
# ==========================
sns.set(style="whitegrid")

# 1. Percent identity histogram
plt.figure(figsize=(8,5))
sns.histplot(df["pident"], bins=50, kde=True, color="skyblue")
plt.title("Percent Identity Distribution")
plt.xlabel("Percent Identity (%)")
plt.ylabel("Number of Hits")
plt.savefig(f"figures/blast/{output_prefix}_pident_hist.png", dpi=300)
plt.close()

# 2. Coverage histogram
plt.figure(figsize=(8,5))
sns.histplot(df["min_cov"], bins=50, kde=True, color="lightgreen")
plt.title("Minimum Coverage Distribution (query vs subject)")
plt.xlabel("Coverage (%)")
plt.ylabel("Number of Hits")
plt.savefig(f"figures/blast/{output_prefix}_coverage_hist.png", dpi=300)
plt.close()

# 3. Identity vs Coverage scatter
plt.figure(figsize=(8,6))
sns.scatterplot(x="pident", y="min_cov", data=df, alpha=0.3)
plt.title("Percent Identity vs Minimum Coverage")
plt.xlabel("Percent Identity (%)")
plt.ylabel("Minimum Coverage (%)")
plt.savefig(f"figures/blast/{output_prefix}_pident_vs_cov.png", dpi=300)
plt.close()

print("Plots generated successfully!")
