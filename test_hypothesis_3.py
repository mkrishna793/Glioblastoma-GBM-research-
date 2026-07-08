import pandas as pd
import numpy as np
import os
import scipy.stats as stats

data_dir = r"D:\new data of the GBM research"
matched_barcodes_path = r"D:\research of the GBM\key_variants_with_barcodes.txt"
out_path = r"D:\research of the GBM\hypothesis_3_results.txt"

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

# Define cohorts
all_rna_patients = set(df_rna_meta['patient_id'].unique())
atrx_mut_rna_patients = sorted(atrx_mut_patients.intersection(all_rna_patients))
atrx_wt_rna_patients = sorted(all_rna_patients.difference(atrx_mut_patients))

# Define gene markers for correlation
mac_markers = ['CD163', 'CD14', 'CSF1R']
inv_markers = ['CHI3L1', 'MMP2', 'MMP9', 'FN1']

# Calculate log2(TPM + 1)
df_log = np.log2(df_rna + 1)

# Perform correlation analysis at recurrence (R1)
tp_target = 'R1'
tp_meta = df_rna_meta[df_rna_meta['timepoint'] == tp_target]

mut_samples = tp_meta[tp_meta['patient_id'].isin(atrx_mut_rna_patients)]['sample'].tolist()
wt_samples = tp_meta[tp_meta['patient_id'].isin(atrx_wt_rna_patients)]['sample'].tolist()

print(f"Timepoint {tp_target} samples: ATRX-Mut = {len(mut_samples)}, ATRX-WT = {len(wt_samples)}")

with open(out_path, "w") as f:
    f.write("="*80 + "\n")
    f.write("HYPOTHESIS 3 TESTING: GENOTYPE-SPECIFIC IMMUNE RECIPROCAL FEEDBACK (ATRX FORK)\n")
    f.write("="*80 + "\n\n")
    
    f.write(f"Analyzed timepoint: {tp_target}\n")
    f.write(f"ATRX-Mutant cohort size: {len(mut_samples)} samples ({len(atrx_mut_rna_patients)} patients)\n")
    f.write(f"ATRX-Wildtype cohort size: {len(wt_samples)} samples ({len(atrx_wt_rna_patients)} patients)\n\n")
    
    f.write("Spearman Correlation between Macrophage Markers and Invasion/Matrix Markers:\n\n")
    
    for mac in mac_markers:
        for inv in inv_markers:
            if mac in df_log.index and inv in df_log.index:
                # ATRX-WT correlation
                wt_vals_mac = df_log.loc[mac, wt_samples].values
                wt_vals_inv = df_log.loc[inv, wt_samples].values
                wt_corr, wt_p = stats.spearmanr(wt_vals_mac, wt_vals_inv)
                
                # ATRX-Mut correlation
                mut_vals_mac = df_log.loc[mac, mut_samples].values
                mut_vals_inv = df_log.loc[inv, mut_samples].values
                
                if len(mut_samples) > 2:
                    mut_corr, mut_p = stats.spearmanr(mut_vals_mac, mut_vals_inv)
                else:
                    mut_corr, mut_p = np.nan, np.nan
                
                f.write(f"Pair: {mac} vs. {inv}\n")
                f.write(f"  ATRX-WT: r = {wt_corr:>6.3f} (p = {wt_p:.2e})\n")
                f.write(f"  ATRX-Mut: r = {mut_corr:>6.3f} (p = {mut_p:.4e} if calculated)\n\n")

print(f"Results saved to {out_path}")
# Print summary to stdout
print("\nCorrelation pairs comparison (ATRX-WT vs ATRX-Mut):")
for mac in ['CD163']:
    for inv in ['CHI3L1', 'MMP2']:
        wt_corr, _ = stats.spearmanr(df_log.loc[mac, wt_samples], df_log.loc[inv, wt_samples])
        mut_corr, _ = stats.spearmanr(df_log.loc[mac, mut_samples], df_log.loc[inv, mut_samples])
        print(f"  {mac} vs {inv} | WT r = {wt_corr:.3f} | Mut r = {mut_corr:.3f}")
