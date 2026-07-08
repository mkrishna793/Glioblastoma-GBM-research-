import pandas as pd
import os

data_dir = r"D:\research of the GBM"

# 1. Load WES Clinical (286 patients)
df_wes_clin = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286_clinical.20200506.txt"), sep="\t")
df_wes_clin.columns = df_wes_clin.columns.str.strip()

print("WES Cohort (286 patients) Summary:")
print("Histology counts:")
print(df_wes_clin['Histology'].value_counts())
print("\nGrade counts:")
print(df_wes_clin['Grade'].value_counts())
print("\nGender counts:")
print(df_wes_clin['Gender'].value_counts())
print("\nTMZ Treatment counts:")
print(df_wes_clin['Chemo_status (TMZ treated=1;un-treated=0)'].value_counts(dropna=False))
print("\nAge Stats:")
print(df_wes_clin['Age'].describe())

# 2. Load IDH-A Clinical (35 patients)
df_idh_clin = pd.read_excel(os.path.join(data_dir, "CGGA_IDH_A_Clinical_Information.xlsx"))
df_idh_clin.columns = df_idh_clin.columns.str.strip()

print("\n" + "="*40)
print("IDH-A Cohort (35 patients) Summary:")
print("Histology counts:")
print(df_idh_clin['Histology'].value_counts())
print("\nGrade counts:")
print(df_idh_clin['Grade'].value_counts())
print("\nGender counts:")
print(df_idh_clin['Gender'].value_counts())
print("\nTMZ Treatment counts:")
print(df_idh_clin['Chemo_status (TMZ treated=1;un-treated=0)'].value_counts(dropna=False))
print("\nAge Stats:")
print(df_idh_clin['Age'].describe())
