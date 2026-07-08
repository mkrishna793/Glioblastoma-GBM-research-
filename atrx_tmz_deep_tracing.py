import pandas as pd
import numpy as np
import os

data_dir = r"D:\research of the GBM"
output_file = os.path.join(data_dir, "atrx_tmz_deep_tracing_report.txt")

# Load datasets
df_mut = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286.20200506.txt"), sep="\t", index_col=0, low_memory=False)
df_mut.columns = df_mut.columns.str.strip()
df_mut_bin = df_mut.transpose().notna().astype(int)
df_mut_bin.index = df_mut_bin.index.astype(str).str.strip()

df_wes_clin = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286_clinical.20200506.txt"), sep="\t")
df_wes_clin.columns = df_wes_clin.columns.str.strip()
df_wes_clin.set_index('CGGA_ID', inplace=True)
df_wes_clin.index = df_wes_clin.index.astype(str).str.strip()

df_rna = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_RNAseq_RSEM_20250915.txt"), sep="\t", index_col=0)
df_rna.columns = df_rna.columns.str.strip()

df_meth = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Methylation_EPIC_Array_20250915.txt"), sep="\t", index_col=0)
df_meth.columns = df_meth.columns.str.strip()

df_prot = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Proteomics_MS_Abundance_20250915.txt"), sep="\t")
df_prot.columns = df_prot.columns.str.strip()

df_phos = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Phosphoproteomics_MS_Abundance_20250915.txt"), sep="\t")
df_phos.columns = df_phos.columns.str.strip()

# Overlap patients
overlap = ['CGGA_P137', 'CGGA_P151', 'CGGA_P153', 'CGGA_P173', 'CGGA_P179', 'CGGA_P271', 'CGGA_P411', 'CGGA_P568', 'CGGA_P633', 'CGGA_P87']

# Grouping
g_wt_treated = ['CGGA_P153', 'CGGA_P179', 'CGGA_P411', 'CGGA_P633']  # ATRX WT, TMZ Treated
g_mut_untreated = ['CGGA_P137', 'CGGA_P173', 'CGGA_P271']          # ATRX Mut, TMZ Untreated
g_mut_treated = ['CGGA_P151', 'CGGA_P568', 'CGGA_P87']              # ATRX Mut, TMZ Treated

report = []
report.append("=" * 80)
report.append("DEEP MOLECULAR TRACING: ATRX MUTATION & TMZ TREATMENT STATUS")
report.append("=" * 80)
report.append(f"Analyzing combinations of ATRX mutation status and TMZ chemotherapy status:")
report.append(f"  Group 1 (ATRX WT, TMZ Treated): N=4 patients ({', '.join(g_wt_treated)})")
report.append(f"  Group 2 (ATRX Mut, TMZ Treated): N=3 patients ({', '.join(g_mut_treated)})")
report.append(f"  Group 3 (ATRX Mut, TMZ Untreated): N=3 patients ({', '.join(g_mut_untreated)})")

# ==============================================================================
# 1. WHAT CHANGES DOES TREATMENT CAUSE IN ATRX MUTANTS? (Group 2 vs Group 3)
# ==============================================================================
report.append("\n" + "=" * 80)
report.append("PART 1: MOLECULAR CHANGES INDUCED BY TMZ TREATMENT IN ATRX MUTANT PATIENTS")
report.append("(Comparing Group 2: ATRX Mut + TMZ Treated vs Group 3: ATRX Mut + TMZ Untreated)")
report.append("=" * 80)

# A. Epigenomics (DNA Methylation)
mean_g2_meth = df_meth[g_mut_treated].mean(axis=1)
mean_g3_meth = df_meth[g_mut_untreated].mean(axis=1)
delta_meth_t = mean_g2_meth - mean_g3_meth

report.append("\n--- Epigenomics (DNA Methylation changes caused by treatment) ---")
report.append("  Top 10 Hypermethylated Probes in Treated (Beta increases):")
for probe, delta in delta_meth_t.sort_values(ascending=False).head(10).items():
    report.append(f"    {probe}: Delta Beta = {delta:.4f} (WT Mean={mean_g3_meth[probe]:.4f} -> Treated Mean={mean_g2_meth[probe]:.4f})")

report.append("  Top 10 Hypomethylated Probes in Treated (Beta decreases):")
for probe, delta in delta_meth_t.sort_values(ascending=True).head(10).items():
    report.append(f"    {probe}: Delta Beta = {delta:.4f} (WT Mean={mean_g3_meth[probe]:.4f} -> Treated Mean={mean_g2_meth[probe]:.4f})")

# B. Transcriptomics (RNA-Seq)
mean_g2_rna = df_rna[g_mut_treated].mean(axis=1)
mean_g3_rna = df_rna[g_mut_untreated].mean(axis=1)
log2fc_rna_t = np.log2((mean_g2_rna + 1) / (mean_g3_rna + 1))

report.append("\n--- Transcriptomics (RNA Expression changes caused by treatment) ---")
report.append("  Top 10 Upregulated Genes in Treated:")
for g, fc in log2fc_rna_t.sort_values(ascending=False).head(10).items():
    report.append(f"    {g}: Log2FC = {fc:.2f}")

report.append("  Top 10 Downregulated Genes in Treated:")
for g, fc in log2fc_rna_t.sort_values(ascending=True).head(10).items():
    report.append(f"    {g}: Log2FC = {fc:.2f}")

# Target Genes Check
target_genes = ['MGMT', 'STAT3', 'EGFR', 'GRIA2', 'NLGN3', 'CD44', 'OLIG2', 'SOX2', 'CHI3L1', 'VIM']
report.append("\n  Expression changes in key markers:")
for tg in target_genes:
    if tg in log2fc_rna_t.index:
        report.append(f"    {tg}: Log2FC = {log2fc_rna_t[tg]:.2f} (Untreated Mean={mean_g3_rna[tg]:.2f} -> Treated Mean={mean_g2_rna[tg]:.2f})")

# C. Proteomics (Protein Abundance)
report.append("\n--- Proteomics (Protein Abundance changes caused by treatment) ---")
for tp in target_genes:
    tp_rows = df_prot[df_prot['Genes'] == tp]
    if not tp_rows.empty:
        row = tp_rows.iloc[0]
        g2_vals = pd.to_numeric(row[g_mut_treated], errors='coerce').mean()
        g3_vals = pd.to_numeric(row[g_mut_untreated], errors='coerce').mean()
        report.append(f"    {tp} Protein: Untreated Mean={g3_vals:.2f} -> Treated Mean={g2_vals:.2f} (Delta={g2_vals-g3_vals:.2f})")

# D. Phosphoproteomics
report.append("\n--- Phosphoproteomics (Active Signaling changes caused by treatment) ---")
target_phos = ['STAT3_Y705', 'STAT3_S727', 'EGFR_S991', 'EGFR_T693', 'AKT1_S129', 'MAPK1_Y187']
for tp in target_phos:
    tp_rows = df_phos[df_phos['PTM.CollapseKey'] == tp]
    if not tp_rows.empty:
        row = tp_rows.iloc[0]
        g2_vals = pd.to_numeric(row[g_mut_treated], errors='coerce').mean()
        g3_vals = pd.to_numeric(row[g_mut_untreated], errors='coerce').mean()
        report.append(f"    {tp} Phospho: Untreated Mean={g3_vals:.2f} -> Treated Mean={g2_vals:.2f} (Delta={g2_vals-g3_vals:.2f})")


# ==============================================================================
# 2. HOW DOES ATRX MUTATION IMPACT TREATED PATIENTS? (Group 2 vs Group 1)
# ==============================================================================
report.append("\n\n" + "=" * 80)
report.append("PART 2: HOW ATRX MUTATION SHAPES THERAPY RESPONSE")
report.append("(Comparing Group 2: ATRX Mut + TMZ Treated vs Group 1: ATRX WT + TMZ Treated)")
report.append("=" * 80)

# A. Epigenomics (DNA Methylation)
mean_g1_meth = df_meth[g_wt_treated].mean(axis=1)
delta_meth_m = mean_g2_meth - mean_g1_meth

report.append("\n--- Epigenomics (DNA Methylation differences in treated patients) ---")
report.append("  Top 10 Hypermethylated Probes in ATRX-Mutant:")
for probe, delta in delta_meth_m.sort_values(ascending=False).head(10).items():
    report.append(f"    {probe}: Delta Beta = {delta:.4f} (WT Mean={mean_g1_meth[probe]:.4f} -> Mut Mean={mean_g2_meth[probe]:.4f})")

report.append("  Top 10 Hypomethylated Probes in ATRX-Mutant:")
for probe, delta in delta_meth_m.sort_values(ascending=True).head(10).items():
    report.append(f"    {probe}: Delta Beta = {delta:.4f} (WT Mean={mean_g1_meth[probe]:.4f} -> Mut Mean={mean_g2_meth[probe]:.4f})")

# B. Transcriptomics (RNA-Seq)
mean_g1_rna = df_rna[g_wt_treated].mean(axis=1)
log2fc_rna_m = np.log2((mean_g2_rna + 1) / (mean_g1_rna + 1))

report.append("\n--- Transcriptomics (RNA Expression differences in treated patients) ---")
report.append("  Top 10 Upregulated Genes in ATRX-Mutant:")
for g, fc in log2fc_rna_m.sort_values(ascending=False).head(10).items():
    report.append(f"    {g}: Log2FC = {fc:.2f}")

report.append("  Top 10 Downregulated Genes in ATRX-Mutant:")
for g, fc in log2fc_rna_m.sort_values(ascending=True).head(10).items():
    report.append(f"    {g}: Log2FC = {fc:.2f}")

report.append("\n  Expression differences in key markers:")
for tg in target_genes:
    if tg in log2fc_rna_m.index:
        report.append(f"    {tg}: Log2FC = {log2fc_rna_m[tg]:.2f} (WT Mean={mean_g1_rna[tg]:.2f} -> Mut Mean={mean_g2_rna[tg]:.2f})")

# C. Proteomics
report.append("\n--- Proteomics (Protein Abundance differences in treated patients) ---")
for tp in target_genes:
    tp_rows = df_prot[df_prot['Genes'] == tp]
    if not tp_rows.empty:
        row = tp_rows.iloc[0]
        g1_vals = pd.to_numeric(row[g_wt_treated], errors='coerce').mean()
        g2_vals = pd.to_numeric(row[g_mut_treated], errors='coerce').mean()
        report.append(f"    {tp} Protein: WT Mean={g1_vals:.2f} -> Mut Mean={g2_vals:.2f} (Delta={g2_vals-g1_vals:.2f})")

# D. Phosphoproteomics
report.append("\n--- Phosphoproteomics (Active Signaling differences in treated patients) ---")
for tp in target_phos:
    tp_rows = df_phos[df_phos['PTM.CollapseKey'] == tp]
    if not tp_rows.empty:
        row = tp_rows.iloc[0]
        g1_vals = pd.to_numeric(row[g_wt_treated], errors='coerce').mean()
        g2_vals = pd.to_numeric(row[g_mut_treated], errors='coerce').mean()
        report.append(f"    {tp} Phospho: WT Mean={g1_vals:.2f} -> Mut Mean={g2_vals:.2f} (Delta={g2_vals-g1_vals:.2f})")

# Save report
with open(output_file, 'w') as f:
    f.write("\n".join(report))

print("Report saved successfully.")
