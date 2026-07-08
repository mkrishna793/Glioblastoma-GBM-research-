import pandas as pd
import numpy as np
import os
from scipy.stats import mannwhitneyu

data_dir = r"D:\research of the GBM"
output_file = os.path.join(data_dir, "idh_mutant_progression_report.txt")

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
df_merged_clin = df_idh_clin.join(df_wes_clin[['Chemo_status (TMZ treated=1;un-treated=0)', 'OS', 'Censor (alive=0; dead=1)']], how='left')

# Load RNA, Methylation, Proteomics, Phosphoproteomics
df_rna = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_RNAseq_RSEM_20250915.txt"), sep="\t", index_col=0)
df_rna.columns = df_rna.columns.str.strip()

df_meth = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Methylation_EPIC_Array_20250915.txt"), sep="\t", index_col=0)
df_meth.columns = df_meth.columns.str.strip()

df_prot = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Proteomics_MS_Abundance_20250915.txt"), sep="\t")
df_prot.columns = df_prot.columns.str.strip()

df_phos = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Phosphoproteomics_MS_Abundance_20250915.txt"), sep="\t")
df_phos.columns = df_phos.columns.str.strip()

# Get patient groupings
g2_patients = df_merged_clin[df_merged_clin['Grade_2021'] == 'G2'].index.tolist()
g4_patients = df_merged_clin[df_merged_clin['Grade_2021'] == 'G4'].index.tolist()

g2_patients = [p for p in g2_patients if p in df_rna.columns]
g4_patients = [p for p in g4_patients if p in df_rna.columns]

chemo_col = 'Chemo_status (TMZ treated=1;un-treated=0)'
treated_patients = df_merged_clin[df_merged_clin[chemo_col] == 1].index.tolist()
untreated_patients = df_merged_clin[df_merged_clin[chemo_col] == 0].index.tolist()

treated_patients = [p for p in treated_patients if p in df_rna.columns]
untreated_patients = [p for p in untreated_patients if p in df_rna.columns]

report = []
report.append("=" * 80)
report.append("IDH-MUTANT CLINICAL PROGRESSION & THERAPY RESPONSE REPORT (N=35)")
report.append("=" * 80)
report.append(f"Cohort size: 35 patients")
report.append(f"Grade 2 (Low Grade): N={len(g2_patients)}")
report.append(f"Grade 4 (High Grade / Aggressive): N={len(g4_patients)}")
report.append(f"TMZ Treated: N={len(treated_patients)}")
report.append(f"TMZ Untreated: N={len(untreated_patients)}")

# ==============================================================================
# SECTION 1: MOLECULAR SWITCHES OF PROGRESSION: GRADE 2 VS GRADE 4
# ==============================================================================
report.append("\n\n" + "=" * 80)
report.append("SECTION 1: MOLECULAR SWITCHES DRIVING GRADE 2 TO GRADE 4 PROGRESSION")
report.append("=" * 80)

# A. Epigenomics (Methylation)
mean_g2_meth = df_meth[g2_patients].mean(axis=1)
mean_g4_meth = df_meth[g4_patients].mean(axis=1)
delta_meth = mean_g4_meth - mean_g2_meth

report.append("\n--- DNA Methylation Changes in Grade 4 Progression ---")
report.append("  Top 10 Hypermethylated Probes in Grade 4 (Loss of gene access):")
for probe, delta in delta_meth.sort_values(ascending=False).head(10).items():
    report.append(f"    {probe}: Delta Beta = {delta:.4f} (G2 Mean={mean_g2_meth[probe]:.4f} -> G4 Mean={mean_g4_meth[probe]:.4f})")

report.append("  Top 10 Hypomethylated Probes in Grade 4 (Gain of gene access):")
for probe, delta in delta_meth.sort_values(ascending=True).head(10).items():
    report.append(f"    {probe}: Delta Beta = {delta:.4f} (G2 Mean={mean_g2_meth[probe]:.4f} -> G4 Mean={mean_g4_meth[probe]:.4f})")

# Check MGMT methylation
mgmt_probes = ['cg12434587', 'cg12981137', 'cg15765353', 'cg21862320']
report.append("\n  MGMT Promoter Methylation in G2 vs G4:")
for probe in mgmt_probes:
    if probe in df_meth.index:
        report.append(f"    {probe}: G2 Mean={df_meth.loc[probe, g2_patients].mean():.4f}, G4 Mean={df_meth.loc[probe, g4_patients].mean():.4f}")

# B. Transcriptomics (RNA)
mean_g2_rna = df_rna[g2_patients].mean(axis=1)
mean_g4_rna = df_rna[g4_patients].mean(axis=1)
log2fc_prog = np.log2((mean_g4_rna + 1) / (mean_g2_rna + 1))

report.append("\n--- RNA Expression Changes in Grade 4 Progression ---")
report.append("  Top 10 Upregulated Genes in Grade 4:")
for g, fc in log2fc_prog.sort_values(ascending=False).head(10).items():
    report.append(f"    {g}: Log2FC = {fc:.2f}")

report.append("  Top 10 Downregulated Genes in Grade 4:")
for g, fc in log2fc_prog.sort_values(ascending=True).head(10).items():
    report.append(f"    {g}: Log2FC = {fc:.2f}")

target_markers = ['MGMT', 'STAT3', 'EGFR', 'GRIA2', 'NLGN3', 'CD44', 'OLIG2', 'SOX2', 'CHI3L1', 'VIM']
report.append("\n  Key Marker Expression in G2 vs G4:")
for tg in target_markers:
    if tg in log2fc_prog.index:
        report.append(f"    {tg}: Log2FC = {log2fc_prog[tg]:.2f} (G2 Mean={mean_g2_rna[tg]:.2f} -> G4 Mean={mean_g4_rna[tg]:.2f})")

# C. Proteomics
report.append("\n--- Protein Abundance Differences in Grade 4 Progression ---")
prot_patient_cols_g2 = [c for c in g2_patients if c in df_prot.columns]
prot_patient_cols_g4 = [c for c in g4_patients if c in df_prot.columns]

for tp in target_markers:
    tp_rows = df_prot[df_prot['Genes'] == tp]
    if not tp_rows.empty:
        row = tp_rows.iloc[0]
        g2_vals = pd.to_numeric(row[prot_patient_cols_g2], errors='coerce').dropna()
        g4_vals = pd.to_numeric(row[prot_patient_cols_g4], errors='coerce').dropna()
        if len(g2_vals) >= 2 and len(g4_vals) >= 2:
            g2_mean = g2_vals.mean()
            g4_mean = g4_vals.mean()
            try:
                stat, pval = mannwhitneyu(g2_vals, g4_vals, alternative='two-sided')
                sig = "SIGNIFICANT" if pval < 0.05 else "not significant"
                report.append(f"    {tp} Protein: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}, p={pval:.4f} ({sig})")
            except:
                report.append(f"    {tp} Protein: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}")

# D. Phosphoproteomics
report.append("\n--- Phospho-Signaling Activation in Grade 4 Progression ---")
phos_patient_cols_g2 = [c for c in g2_patients if c in df_phos.columns]
phos_patient_cols_g4 = [c for c in g4_patients if c in df_phos.columns]
target_phos = ['STAT3_Y705', 'STAT3_S727', 'EGFR_S991', 'EGFR_T693', 'AKT1_S129', 'MAPK1_Y187']

for tp in target_phos:
    tp_rows = df_phos[df_phos['PTM.CollapseKey'] == tp]
    if not tp_rows.empty:
        row = tp_rows.iloc[0]
        g2_vals = pd.to_numeric(row[phos_patient_cols_g2], errors='coerce').dropna()
        g4_vals = pd.to_numeric(row[phos_patient_cols_g4], errors='coerce').dropna()
        if len(g2_vals) >= 2 and len(g4_vals) >= 2:
            g2_mean = g2_vals.mean()
            g4_mean = g4_vals.mean()
            try:
                stat, pval = mannwhitneyu(g2_vals, g4_vals, alternative='two-sided')
                sig = "SIGNIFICANT" if pval < 0.05 else "not significant"
                report.append(f"    {tp} Phospho: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}, p={pval:.4f} ({sig})")
            except:
                report.append(f"    {tp} Phospho: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}")


# ==============================================================================
# SECTION 2: PROTEOMIC SUBTYPES & THEIR MOLECULAR DRIVERS
# ==============================================================================
report.append("\n\n" + "=" * 80)
report.append("SECTION 2: MOLECULAR DRIVERS OF PROTEOMIC CLUSTERS (NEU, AFM, PPR, IME)")
report.append("=" * 80)

# We will group all 35 patients by their Protein_Cluster
clusters = ['NEU', 'AFM', 'PPR', 'IME']
for clust in clusters:
    clust_patients = df_merged_clin[df_merged_clin['Protein_Cluster'] == clust].index.tolist()
    clust_patients = [p for p in clust_patients if p in df_rna.columns]
    
    report.append(f"\n--- Subtype: {clust} (N={len(clust_patients)}) ---")
    
    # Check marker expression
    mean_clust_rna = df_rna[clust_patients].mean(axis=1)
    report.append(f"  Mean RNA expression in {clust}:")
    for tg in target_markers:
        if tg in mean_clust_rna.index:
            report.append(f"    {tg}: {mean_clust_rna[tg]:.2f}")
            
    # Check protein abundance
    report.append(f"  Mean Protein abundance in {clust}:")
    prot_cols_clust = [c for c in clust_patients if c in df_prot.columns]
    for tp in target_markers:
        tp_rows = df_prot[df_prot['Genes'] == tp]
        if not tp_rows.empty:
            row = tp_rows.iloc[0]
            vals = pd.to_numeric(row[prot_cols_clust], errors='coerce').dropna()
            if len(vals) > 0:
                report.append(f"    {tp} Protein: {vals.mean():.2f}")

# Save report
with open(output_file, 'w') as f:
    f.write("\n".join(report))

print("Progression report saved successfully.")
