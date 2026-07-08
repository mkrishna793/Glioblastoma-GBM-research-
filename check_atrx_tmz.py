import pandas as pd
import os

data_dir = r"D:\research of the GBM"

# Load WES mutations
df_mut = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286.20200506.txt"), sep="\t", index_col=0, low_memory=False)
df_mut.columns = df_mut.columns.str.strip()
df_mut_bin = df_mut.transpose().notna().astype(int)
df_mut_bin.index = df_mut_bin.index.astype(str).str.strip()

# Load WES clinical
df_wes_clin = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286_clinical.20200506.txt"), sep="\t")
df_wes_clin.columns = df_wes_clin.columns.str.strip()
df_wes_clin.set_index('CGGA_ID', inplace=True)
df_wes_clin.index = df_wes_clin.index.astype(str).str.strip()

# Overlapping patients
overlap_all = ['CGGA_P137', 'CGGA_P151', 'CGGA_P153', 'CGGA_P173', 'CGGA_P179', 'CGGA_P271', 'CGGA_P401', 'CGGA_P411', 'CGGA_P568', 'CGGA_P633', 'CGGA_P87']

df_overlap_mut = df_mut_bin.loc[overlap_all]
df_overlap_clin = df_wes_clin.loc[overlap_all]

chemo_col = 'Chemo_status (TMZ treated=1;un-treated=0)'

# Create a summary dataframe
df_summary = pd.DataFrame({
    'ATRX': df_overlap_mut['ATRX'],
    'TMZ': df_overlap_clin[chemo_col]
})

print(df_summary)
print("\nCounts:")
print(df_summary.groupby(['ATRX', 'TMZ']).size())
