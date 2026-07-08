import pandas as pd
import numpy as np
import os

data_dir = r"D:\new data of the GBM research"
matched_barcodes_path = r"D:\research of the GBM\key_variants_with_barcodes.txt"

print("Loading matched variants with VAF data...")
df_vars = pd.read_csv(matched_barcodes_path, sep="\t")

# Standardize aliquot_barcode to extract timepoint (TP, R1, R2)
df_vars['timepoint'] = df_vars['aliquot_barcode'].apply(lambda x: x.split("-")[3] if isinstance(x, str) and len(x.split("-")) > 3 else "UNK")
print(f"Variants by timepoint:\n{df_vars['timepoint'].value_counts()}")

# Pivot the dataframe to get VAF (af) for each patient-variant at TP and R1
# Key: (case_barcode, variant_id, gene_symbol) -> af
df_pivot = df_vars.pivot_table(index=['case_barcode', 'variant_id', 'gene_symbol'], 
                              columns='timepoint', 
                              values='af', 
                              aggfunc='mean').fillna(0).reset_index()

# Keep only columns we need
for col in ['TP', 'R1']:
    if col not in df_pivot.columns:
        df_pivot[col] = 0.0

df_pivot['delta_vaf'] = df_pivot['R1'] - df_pivot['TP']
print(f"Pivoted variants shape: {df_pivot.shape}")

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

# Get matched patients with RNA-seq at both TP and R1
tp_patients = set(df_rna_meta[df_rna_meta['timepoint'] == 'TP']['patient_id'].unique())
r1_patients = set(df_rna_meta[df_rna_meta['timepoint'] == 'R1']['patient_id'].unique())
matched_rna_patients = sorted(tp_patients.intersection(r1_patients))
print(f"Matched RNA-seq patients (TP and R1): {len(matched_rna_patients)}")

# For each matched patient, calculate delta RNA for CHI3L1, CD44, and VIM
delta_rna_data = []
for pid in matched_rna_patients:
    tp_sample = df_rna_meta[(df_rna_meta['patient_id'] == pid) & (df_rna_meta['timepoint'] == 'TP')]['sample'].values[0]
    r1_sample = df_rna_meta[(df_rna_meta['patient_id'] == pid) & (df_rna_meta['timepoint'] == 'R1')]['sample'].values[0]
    
    patient_deltas = {'case_barcode': pid}
    for g in ['CHI3L1', 'CD44', 'VIM']:
        if g in df_rna.index:
            tp_val = df_rna.loc[g, tp_sample]
            r1_val = df_rna.loc[g, r1_sample]
            patient_deltas[f'tp_{g}'] = tp_val
            patient_deltas[f'r1_{g}'] = r1_val
            patient_deltas[f'delta_{g}'] = np.log2((r1_val + 1) / (tp_val + 1))
    delta_rna_data.append(patient_deltas)

df_deltas_rna = pd.DataFrame(delta_rna_data)

# Merge DNA VAF changes with RNA expression changes
# A patient can have multiple mutations, let's merge on case_barcode
df_merged = pd.merge(df_pivot, df_deltas_rna, on='case_barcode', how='inner')
print(f"Merged DNA-RNA dataset shape: {df_merged.shape}")

# Save results
out_path = r"D:\research of the GBM\hypothesis_1_results.txt"
with open(out_path, "w") as f:
    f.write("="*80 + "\n")
    f.write("HYPOTHESIS 1 TESTING: CLONAL SELECTION VS. TRANSCRIPTIONAL PLASTICITY\n")
    f.write("="*80 + "\n\n")
    
    # 1. Overall VAF dynamics
    f.write("1. Overall VAF dynamics for pathogenic driver mutations:\n")
    f.write(df_merged[['gene_symbol', 'TP', 'R1', 'delta_vaf']].groupby('gene_symbol').describe().to_string() + "\n\n")
    
    # 2. Correlating VAF changes with RNA changes
    f.write("2. Correlation between VAF change (DNA) and Expression change (RNA):\n")
    for gene in ['IDH1', 'ATRX', 'TP53']:
        sub = df_merged[df_merged['gene_symbol'] == gene]
        f.write(f"\nGene: {gene} (n = {len(sub)} mutation-patient events)\n")
        if not sub.empty:
            for rna_gene in ['CHI3L1', 'CD44', 'VIM']:
                corr = sub['delta_vaf'].corr(sub[f'delta_{rna_gene}'])
                f.write(f"  Corr(delta_vaf, delta_{rna_gene} expression): r = {corr:.4f}\n")
        else:
            f.write("  No overlapping events found.\n")
            
    # 3. Categorizing Selection vs. Plasticity events
    # We define Selection as a significant shift in driver mutation VAF (|delta_vaf| >= 0.10)
    # We define Plasticity as a significant shift in expression (|delta_RNA| >= 1.0) with stable VAF (|delta_vaf| < 0.10)
    f.write("\n3. Classification of resistance events (Selection vs. Plasticity):\n")
    for rna_gene in ['CHI3L1', 'CD44']:
        f.write(f"\nTarget Resistance Gene: {rna_gene}\n")
        # Find cases where the target gene goes up significantly (delta >= 1.0)
        up_cases = df_merged[df_merged[f'delta_{rna_gene}'] >= 1.0]
        f.write(f"  Patients with significant expression increase (log2FC >= 1.0): {len(up_cases['case_barcode'].unique())}\n")
        if not up_cases.empty:
            selection_cases = up_cases[up_cases['delta_vaf'].abs() >= 0.10]
            plasticity_cases = up_cases[up_cases['delta_vaf'].abs() < 0.10]
            
            f.write(f"    - Clonal Selection cases (|delta_vaf| >= 0.10): {len(selection_cases['case_barcode'].unique())}\n")
            f.write(f"    - Pure Plasticity cases (|delta_vaf| < 0.10): {len(plasticity_cases['case_barcode'].unique())}\n")
            if len(selection_cases) > 0:
                f.write("      Selection patient details:\n")
                f.write(selection_cases[['case_barcode', 'gene_symbol', 'delta_vaf', f'delta_{rna_gene}']].to_string() + "\n")
        else:
            f.write("    - No patients with significant expression increase.\n")

print(f"Results saved to {out_path}")
# Print high level statistics to stdout
print("\nClassification of CHI3L1 (YKL-40) expression increase:")
up_cases = df_merged[df_merged['delta_CHI3L1'] >= 1.0]
print(f"  Total events with CHI3L1 increase >= 1.0: {len(up_cases)}")
if not up_cases.empty:
    selection = up_cases[up_cases['delta_vaf'].abs() >= 0.10]
    plasticity = up_cases[up_cases['delta_vaf'].abs() < 0.10]
    print(f"    Clonal Selection: {len(selection)}")
    print(f"    Plasticity: {len(plasticity)}")
