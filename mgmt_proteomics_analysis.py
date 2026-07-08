import pandas as pd
import numpy as np
import os
from scipy.stats import pearsonr

# Define paths
data_dir = r"D:\research of the GBM"
meth_file = os.path.join(data_dir, "CGGA_IDH_A_Methylation_EPIC_Array_20250915.txt")
prot_file = os.path.join(data_dir, "CGGA_IDH_A_Proteomics_MS_Abundance_20250915.txt")
phos_file = os.path.join(data_dir, "CGGA_IDH_A_Phosphoproteomics_MS_Abundance_20250915.txt")
output_report = os.path.join(data_dir, "mgmt_proteomics_report.txt")

print("Loading Methylation data...")
df_meth = pd.read_csv(meth_file, sep="\t", index_col=0)
df_meth.columns = df_meth.columns.str.strip()

print("Loading Proteomics data...")
df_prot = pd.read_csv(prot_file, sep="\t")
df_prot.columns = df_prot.columns.str.strip()

print("Loading Phosphoproteomics data...")
df_phos = pd.read_csv(phos_file, sep="\t")
df_phos.columns = df_phos.columns.str.strip()

# Check for MGMT in Proteomics
# The columns are: 'Protein.Group', 'Genes', followed by patient IDs
print("Unique genes in Proteomics (first 10):", df_prot['Genes'].head(10).tolist())
mgmt_prot_rows = df_prot[df_prot['Genes'] == 'MGMT']
print(f"Found MGMT in Proteomics: {len(mgmt_prot_rows)} rows")

report = []
report.append("=========================================")
report.append("MGMT METHYLATION VS. PROTEOMICS & PHOSPHOPROTEOMICS ANALYSIS")
report.append("=========================================\n")

if not mgmt_prot_rows.empty:
    # We found MGMT! Let's take the first row
    mgmt_row = mgmt_prot_rows.iloc[0]
    protein_id = mgmt_row['Protein.Group']
    
    # Extract patient columns (excluding 'Protein.Group' and 'Genes')
    patient_cols = [c for c in df_prot.columns if c not in ['Protein.Group', 'Genes']]
    
    # Match patient IDs between Methylation and Proteomics
    overlap_patients = df_meth.columns.intersection(patient_cols)
    print(f"Overlap patients for Meth-Prot: {len(overlap_patients)}")
    
    mgmt_prot_vals = pd.to_numeric(mgmt_row[overlap_patients], errors='coerce')
    
    # Let's correlate MGMT protein abundance with top correlated probes
    target_probes = ['cg15765353', 'cg21862320', 'cg13231624', 'cg17583256', 'cg18447506']
    
    report.append("Correlation of top MGMT methylation probes with MGMT Protein Abundance:")
    report.append("Probe_ID\tPearson_R\tP_Value")
    for probe in target_probes:
        if probe in df_meth.index:
            meth_vals = pd.to_numeric(df_meth.loc[probe, overlap_patients], errors='coerce')
            valid_mask = meth_vals.notna() & mgmt_prot_vals.notna()
            if valid_mask.sum() > 5:
                r_val, p_val = pearsonr(meth_vals[valid_mask], mgmt_prot_vals[valid_mask])
                report.append(f"{probe}\t{r_val:.4f}\t{p_val:.4e}")
            else:
                report.append(f"{probe}\tN/A (insufficient overlap)\tN/A")
else:
    report.append("MGMT protein was not detected in the Proteomics dataset.")
    # Check if there are other key target proteins in the dataset (like EGFR, STAT3, TAZ)
    report.append("\nChecking other target proteins in Proteomics:")
    for target in ['EGFR', 'STAT3', 'TAZ', 'CHI3L1', 'SPP1', 'OLIG2', 'SOX2']:
        target_rows = df_prot[df_prot['Genes'] == target]
        report.append(f"  {target}: {len(target_rows)} rows found")

# Check Phosphoproteomics
report.append("\n\n--- PHOSPHOPROTEOMICS TARGET NODES ---")
# The column in Phosphoproteomics is 'PTM.CollapseKey', which represents Gene_Site (e.g. EIF4EBP1_Y34)
# Let's search for phosphorylation keys of our targets (STAT3, EGFR, TAZ)
print("PTM keys (first 10):", df_phos['PTM.CollapseKey'].head(10).tolist())
for target in ['STAT3', 'EGFR', 'TAZ', 'AKT1', 'MAPK1']:
    target_phos = df_phos[df_phos['PTM.CollapseKey'].str.startswith(target, na=False)]
    report.append(f"\nPhosphorylation sites found for {target}:")
    if not target_phos.empty:
        report.append("  PTM.CollapseKey\tValid_Values_Count")
        patient_cols_phos = [c for c in df_phos.columns if c != 'PTM.CollapseKey']
        for idx, row in target_phos.iterrows():
            valid_count = row[patient_cols_phos].notna().sum()
            report.append(f"  {row['PTM.CollapseKey']}\t{valid_count}")
    else:
        report.append("  None found")

# Save report
report_text = "\n".join(report)
with open(output_report, 'w') as f:
    f.write(report_text)

print("Report successfully written to:", output_report)
