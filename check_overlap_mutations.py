import pandas as pd
import os

data_dir = r"D:\research of the GBM"

# Load WES mutations
df_mut = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286.20200506.txt"), sep="\t", index_col=0, low_memory=False)
df_mut.columns = df_mut.columns.str.strip()
df_mut_bin = df_mut.transpose().notna().astype(int)
df_mut_bin.index = df_mut_bin.index.astype(str).str.strip()

# Overlapping patients
overlap_all = ['CGGA_P137', 'CGGA_P151', 'CGGA_P153', 'CGGA_P173', 'CGGA_P179', 'CGGA_P271', 'CGGA_P401', 'CGGA_P411', 'CGGA_P568', 'CGGA_P633', 'CGGA_P87']

df_overlap_mut = df_mut_bin.loc[overlap_all]

# Find genes that have at least some mutated and some wildtype in these 11 patients
mut_counts = df_overlap_mut.sum(axis=0)
varying_genes = mut_counts[(mut_counts >= 2) & (mut_counts <= 9)].sort_values(ascending=False)

print("Genes with variation in the 11 overlap patients:")
for gene, count in varying_genes.items():
    print(f"  {gene}: Mutated={count}, Wildtype={11-count}")
