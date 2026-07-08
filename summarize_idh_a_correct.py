import pandas as pd
import os

data_dir = r"D:\research of the GBM"

# Load IDH-A Clinical (35 patients)
df_idh_clin = pd.read_excel(os.path.join(data_dir, "CGGA_IDH_A_Clinical_Information.xlsx"))
df_idh_clin.columns = df_idh_clin.columns.str.strip()

print("IDH-A Cohort (35 patients) Summary:")
print("\nCohorts counts:")
print(df_idh_clin['Cohorts'].value_counts())

print("\nProtein_Cluster counts:")
print(df_idh_clin['Protein_Cluster'].value_counts())

print("\nGrade_2021 counts:")
print(df_idh_clin['Grade_2021'].value_counts())

print("\nGenomics_Subtype counts:")
print(df_idh_clin['Genomics_Subtype'].value_counts())

print("\nRNA_Subtype counts:")
print(df_idh_clin['RNA_Subtype'].value_counts())

print("\nMethylation_Class counts:")
print(df_idh_clin['Methylation_Class'].value_counts())

print("\nAge Stats:")
print(df_idh_clin['Age'].describe())

print("\nGender counts:")
print(df_idh_clin['Gender'].value_counts())

# Overlapping assays
print("\nAssays summary:")
for col in ['Proteomics', 'RNAseq', 'Phosphoproteomics', 'DNA Methylation EPIC', 'Single Cell RNAseq']:
    if col in df_idh_clin.columns:
        print(f"  {col}: {df_idh_clin[col].notna().sum()} patients available.")
