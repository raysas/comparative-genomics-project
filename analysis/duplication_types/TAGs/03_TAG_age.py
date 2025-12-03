import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

# Assuming 'args' is an object with 'output_dir' attribute and 'tag_pairs' is a DataFrame

# ks_tag and ks_nontag should be defined as the data for TAG and Non-TAG respectively

# === Density Plot ===
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.kdeplot(ks_tag, label='TAG', fill=True, color="#7A0177", alpha=0.6)
sns.kdeplot(ks_nontag, label='Non-TAG', fill=True, color="#084594", alpha=0.6)
plt.axvline(np.median(ks_tag), color="#7A0177", linestyle='--', label='TAG Median')
plt.axvline(np.median(ks_nontag), color="#084594", linestyle='--', label='Non-TAG Median')
plt.title("Ks Distribution: TAG vs Non-TAG", fontsize=16, fontweight='bold')
plt.xlabel("Ks (synonymous substitutions per site)")
plt.ylabel("Density")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "ks_TAG_vs_nonTAG_density.pdf"))
plt.savefig(os.path.join(args.output_dir, "ks_TAG_vs_nonTAG_density.png"), dpi=300)
plt.close()

# === Orientation Analysis ===
if 'Same_Orientation' in tag_pairs.columns:
    print(f"\n=== Orientation Analysis ===")
    orientation_counts = tag_pairs['Same_Orientation'].value_counts()
    same = orientation_counts.get(True, 0)
    total = orientation_counts.sum()
    print(f"Same orientation: {same} / {total} ({same/total*100:.1f}%)")

    # Chi-square test
    obs = [orientation_counts.get(True, 0), orientation_counts.get(False, 0)]
    chi2, p_chi2, _, _ = chi2_contingency([obs, [total/2, total/2]])
    print(f"Chi-square test: χ²={chi2:.2f}, p-value={p_chi2:.2e}")

    # Bar plot
    plt.figure(figsize=(8, 6))
    tag_pairs['Same_Orientation_str'] = tag_pairs['Same_Orientation'].astype(str)
    sns.countplot(x='Same_Orientation_str', data=tag_pairs, hue='Same_Orientation_str', legend=False,
                  palette={'False': "#084594", 'True': "#7A0177"})
    plt.title(f"TAG Gene Orientation\nχ² test: p = {p_chi2:.2e}", fontsize=14, fontweight='bold')
    plt.xlabel("Orientation")
    plt.ylabel("Number of TAG pairs")
    plt.xticks([0, 1], ['Opposite', 'Same'])
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "TAG_orientation.pdf"))
    plt.savefig(os.path.join(args.output_dir, "TAG_orientation.png"), dpi=300)
    plt.close()