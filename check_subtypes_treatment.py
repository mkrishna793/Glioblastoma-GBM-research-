import pandas as pd
import os

data_dir = r"D:\research of the GBM"

# Load IDH-A clinical (35 patients)
df_idh_clin = pd.read_excel(os.path.join(data_dir, "CGGA_IDH_A_Clinical_Information.xlsx"))
df_idh_clin.columns = df_idh_clin.columns.str.strip()
df_idh_clin.set_index('Cohort_ID', inplace=True)
df_idh_clin.index = df_idh_clin.index.astype(str).str.strip()

# Load WES Clinical (286 patients) to get chemo status
df_wes_clin = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286_clinical.20200506.txt"), sep="\t")
df_wes_clin.columns = df_wes_clin.columns.str.strip()
df_wes_clin.set_index('CGGA_ID', inplace=True)
df_wes_clin.index = df_wes_clin.index.astype(str).str.strip()

# Merge clinical data
df_merged_clin = df_idh_clin.join(df_wes_clin[['Chemo_status (TMZ treated=1;un-treated=0)']], how='left')
df_merged_clin.rename(columns={'Chemo_status (TMZ treated=1;un-treated=0)': 'TMZ'}, inplace=True)

# Count distribution of TMZ (Treated=1.0, Untreated=0.0) across Protein_Cluster
distribution = df_merged_clin.groupby(['Protein_Cluster', 'TMZ'], dropna=False).size().unstack(fill_value=0)
print("Distribution of Subtypes by TMZ Treatment Status:")
print(distribution)
