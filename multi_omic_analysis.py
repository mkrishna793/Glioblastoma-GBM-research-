import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr

# Define paths
data_dir = r"D:\research of the GBM"
clinical_xlsx = os.path.join(data_dir, "CGGA_IDH_A_Clinical_Information.xlsx")
rna_file = os.path.join(data_dir, "CGGA_IDH_A_RNAseq_RSEM_20250915.txt")
meth_file = os.path.join(data_dir, "CGGA_IDH_A_Methylation_EPIC_Array_20250915.txt")
output_report = os.path.join(data_dir, "multi_omic_integration_report.txt")

# WES files
wes_clin_file = os.path.join(data_dir, "CGGA.WEseq_286_clinical.20200506.txt")
wes_mut_file = os.path.join(data_dir, "CGGA.WEseq_286.20200506.txt")

print("Loading datasets...")
df_clin_idh = pd.read_excel(clinical_xlsx)
df_clin_idh.columns = df_clin_idh.columns.str.strip()
df_clin_idh.set_index('Cohort_ID', inplace=True)
df_clin_idh.index = df_clin_idh.index.astype(str).str.strip()

df_rna = pd.read_csv(rna_file, sep="\t", index_col=0)
df_rna.columns = df_rna.columns.str.strip()

# Methylation file is very large (435 MB), let's load it efficiently
df_meth = pd.read_csv(meth_file, sep="\t", index_col=0)
df_meth.columns = df_meth.columns.str.strip()

print("Datasets loaded successfully.")
print(f"IDH-A Clinical: {df_clin_idh.shape}")
print(f"RNA-Seq: {df_rna.shape}")
print(f"Methylation: {df_meth.shape}")

# Find overlapping patients across Clinical, RNA, and Methylation
overlap_patients = df_clin_idh.index.intersection(df_rna.columns).intersection(df_meth.columns)
print(f"Overlapping patients (Clinical + RNA + Methylation): {len(overlap_patients)}")
print(list(overlap_patients))

# Create report buffer
report = []
report.append("=========================================")
report.append("CGGA GLIOBLASTOMA MULTI-OMIC INTEGRATION REPORT")
report.append("=========================================\n")

report.append("--- PART 1: PATIENT MATCHING ACROSS LAYERS ---")
report.append(f"Cohort Size (Clinical): {len(df_clin_idh)}")
report.append(f"Overlapping Patients (Clinical + RNA + Methylation): {len(overlap_patients)}")
report.append(", ".join(list(overlap_patients)))

# Analysis A: MGMT Epigenetic-to-Transcriptomic Translation
report.append("\n\n--- PART 2: MGMT EPIGENETIC-TO-TRANSCRIPTOMIC SILENCING ANALYSIS ---")
if 'MGMT' in df_rna.index:
    mgmt_rna = df_rna.loc['MGMT', overlap_patients]
    
    # We want to find which methylation probes in the EPIC array correlate with MGMT RNA expression.
    # To do this efficiently, we will calculate the correlation between MGMT RNA and EVERY probe in the methylation dataset.
    # And then we will pull the top negatively correlated probes.
    print("Calculating MGMT methylation-expression correlations...")
    correlations = []
    for probe_id, row in df_meth.iterrows():
        # Clean data (drop NaNs)
        meth_vals = row[overlap_patients]
        valid_mask = meth_vals.notna() & mgmt_rna.notna()
        if valid_mask.sum() > 10: # Minimum 10 valid data points
            r_val, p_val = pearsonr(meth_vals[valid_mask], mgmt_rna[valid_mask])
            if not np.isnan(r_val):
                correlations.append((probe_id, r_val, p_val))
                
    df_corr = pd.DataFrame(correlations, columns=['Probe_ID', 'Pearson_R', 'P_Value'])
    df_corr_sorted = df_corr.sort_values(by='Pearson_R', ascending=True) # Strongest negative correlation first
    
    report.append("Top 15 Methylation Probes Negatively Correlated with MGMT Expression (Gene Silencing Switches):")
    report.append("Probe_ID\tPearson_R\tP_Value")
    for idx, row in df_corr_sorted.head(15).iterrows():
        report.append(f"{row['Probe_ID']}\t{row['Pearson_R']:.4f}\t{row['P_Value']:.4e}")
else:
    report.append("MGMT gene not found in RNA-Seq dataset.")

# Analysis B: Treatment-specific Methylation & Expression Shifts
report.append("\n\n--- PART 3: TREATMENT-SPECIFIC METHYLATION & EXPRESSION SHIFTS ---")
# Let's check how Chemotherapy (TMZ) and Radiotherapy affect MGMT expression and its top silencing probe (cg12434587 or whichever is top)
top_probe = df_corr_sorted.iloc[0]['Probe_ID'] if 'df_corr_sorted' in locals() and not df_corr_sorted.empty else None

if top_probe:
    # Match patient IDs and check their Chemotherapy/Radiotherapy status
    # Note: In IDH_A cohort, the clinical spreadsheet has column names: 'OS_day', 'Censor_OS', 'RNA_Subtype', etc.
    # Let's inspect the clinical data to find the chemotherapy or radiotherapy status.
    # In the clinical table columns we saw: 'Cohorts', 'PR_Status', 'Grade_2021', 'RNA_Subtype', 'DNA Methylation EPIC'
    # Wait, does the IDH-A clinical spreadsheet have chemo/radio columns?
    # Let's check the columns: 'Cohort_ID', 'Cohorts', 'PR_Status', 'Grade_2021', 'RNA_Subtype', etc.
    # It does not have Chemo_status directly in the column list, but it has 'PR_Status' (Primary vs Recurrent).
    # However, our WES 286 clinical dataset DOES have 'Chemo_status (TMZ treated=1;un-treated=0)' and 'Radio_status'.
    # So we can run the treatment analysis on the WES 286 dataset!
    report.append(f"Analyzing treatment-specific shifts using WES 286 clinical dataset:")
    df_wes_clin = pd.read_csv(wes_clin_file, sep="\t")
    df_wes_clin.columns = df_wes_clin.columns.str.strip()
    df_wes_clin.set_index('CGGA_ID', inplace=True)
    df_wes_clin.index = df_wes_clin.index.astype(str).str.strip()
    
    # We have WES mutation and clinical data for 286 patients.
    # Let's check if the clinical file has 'MGMTp_meth' (methylated vs un-methylated)
    # We saw in columns: 'MGMTp_methylation_status'
    # Let's compare TMZ Chemotherapy status vs. MGMT methylation status and patient overall survival (OS).
    # This will show the therapeutic dilemma!
    df_wes_clin['OS'] = pd.to_numeric(df_wes_clin['OS'], errors='coerce')
    df_wes_clin['Censor (alive=0; dead=1)'] = pd.to_numeric(df_wes_clin['Censor (alive=0; dead=1)'], errors='coerce')
    df_wes_valid = df_wes_clin.dropna(subset=['OS', 'Censor (alive=0; dead=1)'])
    
    # Group by: Chemo treated (1) vs. untreated (0) and MGMT status
    groups = [
        ("Chemo Treated + MGMT Methylated", df_wes_valid[(df_wes_valid['Chemo_status (TMZ treated=1;un-treated=0)'] == 1) & (df_wes_valid['MGMTp_methylation_status'] == 'methylated')]),
        ("Chemo Treated + MGMT Unmethylated", df_wes_valid[(df_wes_valid['Chemo_status (TMZ treated=1;un-treated=0)'] == 1) & (df_wes_valid['MGMTp_methylation_status'] == 'un-methylated')]),
        ("Chemo Untreated + MGMT Methylated", df_wes_valid[(df_wes_valid['Chemo_status (TMZ treated=1;un-treated=0)'] == 0) & (df_wes_valid['MGMTp_methylation_status'] == 'methylated')]),
        ("Chemo Untreated + MGMT Unmethylated", df_wes_valid[(df_wes_valid['Chemo_status (TMZ treated=1;un-treated=0)'] == 0) & (df_wes_valid['MGMTp_methylation_status'] == 'un-methylated')]),
    ]
    
    report.append("\nSurvival Outcomes Across Chemotherapy & Epigenetic Status (The TMZ Chemosensitivity Loop):")
    report.append("Group\tCount\tMean OS (days)\tDeath Rate (%)")
    for g_name, g_df in groups:
        count = len(g_df)
        mean_os = g_df['OS'].mean() if count > 0 else 0
        death_rate = (g_df['Censor (alive=0; dead=1)'].mean() * 100) if count > 0 else 0
        report.append(f"{g_name}\t{count}\t{mean_os:.1f}\t{death_rate:.1f}%")

# Analysis C: Transcriptional Subtypes vs. Survival (IDH-A cohort)
report.append("\n\n--- PART 4: TRANSCRIPTIONAL SUBTYPES & SURVIVAL OUTCOMES ---")
# RNA_Subtype vs OS_day in IDH-A
df_clin_idh['OS_day'] = pd.to_numeric(df_clin_idh['OS_day'], errors='coerce')
df_clin_idh['Censor_OS'] = pd.to_numeric(df_clin_idh['Censor_OS'], errors='coerce')
df_idh_valid = df_clin_idh.dropna(subset=['OS_day', 'Censor_OS'])

subtypes = df_idh_valid['RNA_Subtype'].unique()
report.append("RNA_Subtype\tCount\tMean OS (days)\tDeath Rate (%)")
for st in subtypes:
    st_df = df_idh_valid[df_idh_valid['RNA_Subtype'] == st]
    count = len(st_df)
    mean_os = st_df['OS_day'].mean() if count > 0 else 0
    death_rate = (st_df['Censor_OS'].mean() * 100) if count > 0 else 0
    report.append(f"{st}\t{count}\t{mean_os:.1f}\t{death_rate:.1f}%")

# Save report
report_text = "\n".join(report)
with open(output_report, 'w') as f:
    f.write(report_text)

print("Report successfully written to:", output_report)
