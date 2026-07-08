import pandas as pd
import os

data_dir = r"D:\research of the GBM"

# Load WES mutations
df_mut = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286.20200506.txt"), sep="\t", index_col=0, low_memory=False)
df_mut.columns = df_mut.columns.str.strip()
df_mut_bin = df_mut.transpose().notna().astype(int)

# Calculate mutation rates
mut_rates = df_mut_bin.mean(axis=0) * 100
top_mutated = mut_rates.sort_values(ascending=False).head(20)

print("Top 20 Mutated Genes in the WES 286 Cohort:")
for gene, rate in top_mutated.items():
    print(f"  {gene}: {rate:.1f}%")
