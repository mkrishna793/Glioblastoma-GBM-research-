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

# Helper to filter patient lists
def align_patients(cohort_list, df_cols):
    return [p for p in cohort_list if p in df_cols]

# Get base patient groupings
g2_all = df_merged_clin[df_merged_clin['Grade_2021'] == 'G2'].index.tolist()
g4_all = df_merged_clin[df_merged_clin['Grade_2021'] == 'G4'].index.tolist()

chemo_col = 'Chemo_status (TMZ treated=1;un-treated=0)'
treated_all = df_merged_clin[df_merged_clin[chemo_col] == 1].index.tolist()
untreated_all = df_merged_clin[df_merged_clin[chemo_col] == 0].index.tolist()

report = []
report.append("=" * 80)
report.append("IDH-MUTANT CLINICAL PROGRESSION & THERAPY RESPONSE REPORT (N=35)")
report.append("=" * 80)
report.append(f"Cohort size: 35 patients")
report.append(f"Grade 2 (Low Grade): N={len(g2_all)}")
report.append(f"Grade 4 (High Grade / Aggressive): N={len(g4_all)}")
report.append(f"TMZ Treated: N={len(treated_all)}")
report.append(f"TMZ Untreated: N={len(untreated_all)}")

# ==============================================================================
# SECTION 1: MOLECULAR SWITCHES OF PROGRESSION: GRADE 2 VS GRADE 4
# ==============================================================================
report.append("\n\n" + "=" * 80)
report.append("SECTION 1: MOLECULAR SWITCHES DRIVING GRADE 2 TO GRADE 4 PROGRESSION")
report.append("=" * 80)

# A. Epigenomics (Methylation)
g2_meth = align_patients(g2_all, df_meth.columns)
g4_meth = align_patients(g4_all, df_meth.columns)

report.append(f"\n--- DNA Methylation Changes in Grade 4 Progression (G2={len(g2_meth)}, G4={len(g4_meth)}) ---")
if len(g2_meth) >= 2 and len(g4_meth) >= 2:
    mean_g2_meth = df_meth[g2_meth].mean(axis=1)
    mean_g4_meth = df_meth[g4_meth].mean(axis=1)
    delta_meth = mean_g4_meth - mean_g2_meth

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
            report.append(f"    {probe}: G2 Mean={df_meth.loc[probe, g2_meth].mean():.4f}, G4 Mean={df_meth.loc[probe, g4_meth].mean():.4f}")
else:
    report.append("  Insufficient overlap for methylation comparison.")

# B. Transcriptomics (RNA)
g2_rna = align_patients(g2_all, df_rna.columns)
g4_rna = align_patients(g4_all, df_rna.columns)

report.append(f"\n--- RNA Expression Changes in Grade 4 Progression (G2={len(g2_rna)}, G4={len(g4_rna)}) ---")
if len(g2_rna) >= 2 and len(g4_rna) >= 2:
    mean_g2_rna = df_rna[g2_rna].mean(axis=1)
    mean_g4_rna = df_rna[g4_rna].mean(axis=1)
    log2fc_prog = np.log2((mean_g4_rna + 1) / (mean_g2_rna + 1))

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
else:
    report.append("  Insufficient overlap for RNA comparison.")

# C. Proteomics
prot_cols = [c for c in df_prot.columns if c not in ['Protein.Group', 'Genes']]
g2_prot = align_patients(g2_all, prot_cols)
g4_prot = align_patients(g4_all, prot_cols)

report.append(f"\n--- Protein Abundance Differences in Grade 4 Progression (G2={len(g2_prot)}, G4={len(g4_prot)}) ---")
if len(g2_prot) >= 2 and len(g4_prot) >= 2:
    for tp in target_markers:
        tp_rows = df_prot[df_prot['Genes'] == tp]
        if not tp_rows.empty:
            row = tp_rows.iloc[0]
            g2_vals = pd.to_numeric(row[g2_prot], errors='coerce').dropna()
            g4_vals = pd.to_numeric(row[g4_prot], errors='coerce').dropna()
            if len(g2_vals) >= 2 and len(g4_vals) >= 2:
                g2_mean = g2_vals.mean()
                g4_mean = g4_vals.mean()
                try:
                    stat, pval = mannwhitneyu(g2_vals, g4_vals, alternative='two-sided')
                    sig = "SIGNIFICANT" if pval < 0.05 else "not significant"
                    report.append(f"    {tp} Protein: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}, p={pval:.4f} ({sig})")
                except Exception as e:
                    report.append(f"    {tp} Protein: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}")
else:
    report.append("  Insufficient overlap for proteomics comparison.")

# D. Phosphoproteomics
phos_cols = [c for c in df_phos.columns if c != 'PTM.CollapseKey']
g2_phos = align_patients(g2_all, phos_cols)
g4_phos = align_patients(g4_all, phos_cols)
target_phos = ['STAT3_Y705', 'STAT3_S727', 'EGFR_S991', 'EGFR_T693', 'AKT1_S129', 'MAPK1_Y187']

report.append(f"\n--- Phospho-Signaling Activation in Grade 4 Progression (G2={len(g2_phos)}, G4={len(g4_phos)}) ---")
if len(g2_phos) >= 2 and len(g4_phos) >= 2:
    for tp in target_phos:
        tp_rows = df_phos[df_phos['PTM.CollapseKey'] == tp]
        if not tp_rows.empty:
            row = tp_rows.iloc[0]
            g2_vals = pd.to_numeric(row[g2_phos], errors='coerce').dropna()
            g4_vals = pd.to_numeric(row[g4_phos], errors='coerce').dropna()
            if len(g2_vals) >= 2 and len(g4_vals) >= 2:
                g2_mean = g2_vals.mean()
                g4_mean = g4_vals.mean()
                try:
                    stat, pval = mannwhitneyu(g2_vals, g4_vals, alternative='two-sided')
                    sig = "SIGNIFICANT" if pval < 0.05 else "not significant"
                    report.append(f"    {tp} Phospho: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}, p={pval:.4f} ({sig})")
                except Exception as e:
                    report.append(f"    {tp} Phospho: G2 Mean={g2_mean:.2f}, G4 Mean={g4_mean:.2f}")
else:
    report.append("  Insufficient overlap for phosphoproteomics comparison.")


# ==============================================================================
# SECTION 2: PROTEOMIC SUBTYPES & THEIR MOLECULAR DRIVERS
# ==============================================================================
report.append("\n\n" + "=" * 80)
report.append("SECTION 2: MOLECULAR DRIVERS OF PROTEOMIC CLUSTERS (NEU, AFM, PPR, IME)")
report.append("=" * 80)

clusters = ['NEU', 'AFM', 'PPR', 'IME']
for clust in clusters:
    clust_patients = df_merged_clin[df_merged_clin['Protein_Cluster'] == clust].index.tolist()
    
    clust_rna = align_patients(clust_patients, df_rna.columns)
    clust_prot = align_patients(clust_patients, prot_cols)
    
    report.append(f"\n--- Subtype: {clust} (N_RNA={len(clust_rna)}, N_Prot={len(clust_prot)}) ---")
    
    if len(clust_rna) > 0:
        mean_clust_rna = df_rna[clust_rna].mean(axis=1)
        report.append(f"  Mean RNA expression in {clust}:")
        for tg in target_markers:
            if tg in mean_clust_rna.index:
                report.append(f"    {tg}: {mean_clust_rna[tg]:.2f}")
            
    if len(clust_prot) > 0:
        report.append(f"  Mean Protein abundance in {clust}:")
        for tp in target_markers:
            tp_rows = df_prot[df_prot['Genes'] == tp]
            if not tp_rows.empty:
                row = tp_rows.iloc[0]
                vals = pd.to_numeric(row[clust_prot], errors='coerce').dropna()
                if len(vals) > 0:
                    report.append(f"    {tp} Protein: {vals.mean():.2f}")

# Save report
with open(output_file, 'w') as f:
    f.write("\n".join(report))

print("Progression report saved successfully.")
