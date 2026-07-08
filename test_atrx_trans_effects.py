import pandas as pd
import numpy as np
import os
from scipy.stats import mannwhitneyu

data_dir = r"D:\new data of the GBM research"
matched_barcodes_path = r"D:\research of the GBM\key_variants_with_barcodes.txt"

print("Loading GLASS RNA-seq TPM data...")
df_rna = pd.read_csv(os.path.join(data_dir, "gene_tpm_matrix_all_samples.tsv"), sep="\t", index_col=0)
df_rna.columns = df_rna.columns.str.strip()

# Standardize RNA columns to patient ID + timepoint
rna_sample_meta = []
for col in df_rna.columns:
    parts = col.replace("-", ".").split(".")
    pid = "-".join(parts[:3]) # Standardize to GLSS-XX-XXXX
    tp = parts[3] if len(parts) > 3 else "UNK"
    rna_sample_meta.append({"sample": col, "patient_id": pid, "timepoint": tp})
df_rna_meta = pd.DataFrame(rna_sample_meta)

print("Loading final matched variants with patient barcodes...")
df_vars = pd.read_csv(matched_barcodes_path, sep="\t")

# Find ATRX mutant patient IDs
atrx_mut_patients = set(df_vars[df_vars['gene_symbol'] == 'ATRX']['case_barcode'].unique())
print(f"ATRX-mutant patients in variants: {len(atrx_mut_patients)} ({atrx_mut_patients})")

# Define cohorts for primary (TP) and recurrence (R1)
all_rna_patients = set(df_rna_meta['patient_id'].unique())
atrx_mut_rna_patients = atrx_mut_patients.intersection(all_rna_patients)
atrx_wt_rna_patients = all_rna_patients.difference(atrx_mut_patients)

print(f"RNA-seq Cohort Breakdown:")
print(f"  ATRX-Mutant Patients: {len(atrx_mut_rna_patients)}")
print(f"  ATRX-Wildtype Patients: {len(atrx_wt_rna_patients)}")

# Analyze targets: CHI3L1, GRIA2, NLGN3, CD44 at TP and R1
targets = ['CHI3L1', 'GRIA2', 'NLGN3', 'CD44']

for tp in ['TP', 'R1']:
    print(f"\n==================== TIMEPOINT: {tp} ====================")
    # Get samples for this timepoint
    tp_meta = df_rna_meta[df_rna_meta['timepoint'] == tp]
    
    mut_samples = tp_meta[tp_meta['patient_id'].isin(atrx_mut_rna_patients)]['sample'].tolist()
    wt_samples = tp_meta[tp_meta['patient_id'].isin(atrx_wt_rna_patients)]['sample'].tolist()
    
    print(f"Samples: ATRX-Mut = {len(mut_samples)}, ATRX-WT = {len(wt_samples)}")
    
    if len(mut_samples) > 0 and len(wt_samples) > 0:
        for g in targets:
            if g in df_rna.index:
                mut_vals = df_rna.loc[g, mut_samples].values
                wt_vals = df_rna.loc[g, wt_samples].values
                
                mut_median = np.median(mut_vals)
                wt_median = np.median(wt_vals)
                
                # Perform Mann-Whitney U test (trans-regulatory association)
                stat, p_val = mannwhitneyu(mut_vals, wt_vals, alternative='two-sided')
                
                log2_fc = np.log2((mut_median + 1) / (wt_median + 1))
                
                print(f"  {g:<8} | WT Med: {wt_median:>7.2f} | Mut Med: {mut_median:>7.2f} | Log2FC: {log2_fc:>5.2f} | MW-p: {p_val:.4f}")
            else:
                print(f"  {g:<8} | Not found in RNA-seq matrix")
