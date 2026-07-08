import pandas as pd
import numpy as np
import os
from scipy.stats import mannwhitneyu

data_dir = r"D:\research of the GBM"
output_file = os.path.join(data_dir, "complete_variant_tracing_report.txt")

print("=" * 60)
print("LOADING ALL CGGA DATASETS")
print("=" * 60)

# 1. WES Mutations (286 patients, 12574 genes)
print("Loading WES mutations...")
df_mut = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286.20200506.txt"), sep="\t", index_col=0, low_memory=False)
df_mut.columns = df_mut.columns.str.strip()
df_mut_bin = df_mut.transpose().notna().astype(int)
df_mut_bin.index = df_mut_bin.index.astype(str).str.strip()

# 2. WES Clinical (286 patients with treatment info)
print("Loading WES clinical...")
df_wes_clin = pd.read_csv(os.path.join(data_dir, "CGGA.WEseq_286_clinical.20200506.txt"), sep="\t")
df_wes_clin.columns = df_wes_clin.columns.str.strip()
df_wes_clin.set_index('CGGA_ID', inplace=True)
df_wes_clin.index = df_wes_clin.index.astype(str).str.strip()
df_wes_clin['OS'] = pd.to_numeric(df_wes_clin['OS'], errors='coerce')
df_wes_clin['Censor (alive=0; dead=1)'] = pd.to_numeric(df_wes_clin['Censor (alive=0; dead=1)'], errors='coerce')

# 3. IDH-A Clinical (35 patients)
print("Loading IDH-A clinical...")
df_idh_clin = pd.read_excel(os.path.join(data_dir, "CGGA_IDH_A_Clinical_Information.xlsx"))
df_idh_clin.columns = df_idh_clin.columns.str.strip()
df_idh_clin.set_index('Cohort_ID', inplace=True)
df_idh_clin.index = df_idh_clin.index.astype(str).str.strip()

# 4. RNA-Seq (genes x patients)
print("Loading RNA-Seq...")
df_rna = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_RNAseq_RSEM_20250915.txt"), sep="\t", index_col=0)
df_rna.columns = df_rna.columns.str.strip()

# 5. Methylation EPIC Array (probes x patients)
print("Loading Methylation...")
df_meth = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Methylation_EPIC_Array_20250915.txt"), sep="\t", index_col=0)
df_meth.columns = df_meth.columns.str.strip()

# 6. Proteomics
print("Loading Proteomics...")
df_prot = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Proteomics_MS_Abundance_20250915.txt"), sep="\t")
df_prot.columns = df_prot.columns.str.strip()

# 7. Phosphoproteomics
print("Loading Phosphoproteomics...")
df_phos = pd.read_csv(os.path.join(data_dir, "CGGA_IDH_A_Phosphoproteomics_MS_Abundance_20250915.txt"), sep="\t")
df_phos.columns = df_phos.columns.str.strip()

# Normal Brain RNA
print("Loading Normal Brain RNA...")
df_normal = pd.read_csv(os.path.join(data_dir, "CGGA.normal_20.Read_Counts-genes.20230104.txt"), sep="\t", index_col=0)
total_counts = df_normal.sum(axis=0)
df_normal_cpm = (df_normal / total_counts) * 1e6

print("\nAll datasets loaded.")

# Find overlapping patients across WES + IDH-A layers
wes_patients = set(df_mut_bin.index)
idh_patients = set(df_idh_clin.index)
rna_patients = set(df_rna.columns)
meth_patients = set(df_meth.columns)
prot_patient_cols = [c for c in df_prot.columns if c not in ['Protein.Group', 'Genes']]
phos_patient_cols = [c for c in df_phos.columns if c != 'PTM.CollapseKey']

overlap_all = wes_patients & idh_patients & rna_patients & meth_patients & set(prot_patient_cols)
overlap_wes_idh = wes_patients & idh_patients
print(f"WES-IDH overlap: {len(overlap_wes_idh)} patients")
print(f"Full overlap (WES+RNA+Meth+Prot): {len(overlap_all)} patients")
print(f"Full overlap patients: {sorted(overlap_all)}")

# =========================================================
# BUILD THE REPORT
# =========================================================
report = []
report.append("=" * 80)
report.append("COMPLETE VARIANT-TO-PROTEOME TRACING REPORT")
report.append("Tracing each mutated gene through: Genomics -> Epigenomics -> Transcriptomics -> Proteomics")
report.append("=" * 80)
report.append(f"\nWES Cohort: {len(df_mut_bin)} patients")
report.append(f"IDH-A Cohort: {len(df_idh_clin)} patients")
report.append(f"Full Multi-Omic Overlap: {len(overlap_all)} patients")
report.append(f"Overlap Patient IDs: {', '.join(sorted(overlap_all))}")

# =========================================================
# TARGET GENES TO TRACE
# =========================================================
target_genes = ['IDH1', 'TP53', 'ATRX', 'PTEN', 'EGFR', 'CIC', 'NF1', 'NOTCH1', 'PIK3CA']

for gene in target_genes:
    report.append("\n" + "=" * 80)
    report.append(f"TRACING GENE: {gene}")
    report.append("=" * 80)
    
    if gene not in df_mut_bin.columns:
        report.append(f"  {gene} not found in WES mutation matrix. Skipping.")
        continue
    
    # ---- LAYER 1: GENOMICS (Mutation Status) ----
    report.append(f"\n--- LAYER 1: GENOMICS (WES 286 Cohort) ---")
    mutated_patients = df_mut_bin.index[df_mut_bin[gene] == 1].tolist()
    wildtype_patients = df_mut_bin.index[df_mut_bin[gene] == 0].tolist()
    report.append(f"  Mutated: {len(mutated_patients)} patients")
    report.append(f"  Wildtype: {len(wildtype_patients)} patients")
    report.append(f"  Mutation Rate: {len(mutated_patients)/len(df_mut_bin)*100:.1f}%")
    
    # Survival by mutation status
    for group_name, group_ids in [("Mutated", mutated_patients), ("Wildtype", wildtype_patients)]:
        group_clin = df_wes_clin.loc[df_wes_clin.index.isin(group_ids)].dropna(subset=['OS'])
        if len(group_clin) > 0:
            mean_os = group_clin['OS'].mean()
            death_rate = group_clin['Censor (alive=0; dead=1)'].mean() * 100
            report.append(f"  {group_name}: Mean OS = {mean_os:.0f} days, Death Rate = {death_rate:.1f}%")
    
    # Treatment stratification
    report.append(f"\n  Treatment Stratification for {gene}:")
    chemo_col = 'Chemo_status (TMZ treated=1;un-treated=0)'
    for mut_status, mut_label in [(1, "Mutated"), (0, "Wildtype")]:
        for chemo_status, chemo_label in [(1, "TMZ Treated"), (0, "TMZ Untreated")]:
            subset_ids = df_mut_bin.index[df_mut_bin[gene] == mut_status].tolist()
            subset_clin = df_wes_clin.loc[df_wes_clin.index.isin(subset_ids)]
            subset_clin = subset_clin[subset_clin[chemo_col] == chemo_status].dropna(subset=['OS'])
            if len(subset_clin) >= 3:
                mean_os = subset_clin['OS'].mean()
                death_rate = subset_clin['Censor (alive=0; dead=1)'].mean() * 100
                report.append(f"    {mut_label} + {chemo_label}: N={len(subset_clin)}, Mean OS={mean_os:.0f} days, Death Rate={death_rate:.1f}%")
    
    # ---- LAYER 2: EPIGENOMICS (Methylation) ----
    report.append(f"\n--- LAYER 2: EPIGENOMICS (DNA Methylation) ---")
    # For patients in the overlap, compare methylation between mutated and wildtype
    mut_overlap = [p for p in mutated_patients if p in meth_patients]
    wt_overlap = [p for p in wildtype_patients if p in meth_patients]
    
    if len(mut_overlap) >= 2 and len(wt_overlap) >= 2:
        report.append(f"  Comparing methylation: {len(mut_overlap)} mutated vs {len(wt_overlap)} wildtype patients")
        
        # Calculate mean methylation for each group and find differentially methylated probes
        mean_mut_meth = df_meth[mut_overlap].mean(axis=1)
        mean_wt_meth = df_meth[wt_overlap].mean(axis=1)
        delta_meth = mean_mut_meth - mean_wt_meth
        
        # Top hypermethylated probes in mutated group
        top_hyper = delta_meth.sort_values(ascending=False).head(10)
        report.append(f"  Top 10 Hypermethylated Probes (higher in {gene}-mutated):")
        for probe, delta in top_hyper.items():
            report.append(f"    {probe}: Delta Beta = {delta:.4f}")
        
        # Top hypomethylated probes in mutated group
        top_hypo = delta_meth.sort_values(ascending=True).head(10)
        report.append(f"  Top 10 Hypomethylated Probes (lower in {gene}-mutated):")
        for probe, delta in top_hypo.items():
            report.append(f"    {probe}: Delta Beta = {delta:.4f}")
        
        # Check MGMT promoter region probes
        report.append(f"  MGMT-related Methylation in {gene}-mutated vs wildtype:")
        mgmt_probes = ['cg12434587', 'cg12981137', 'cg15765353', 'cg21862320']
        for probe in mgmt_probes:
            if probe in df_meth.index:
                mut_val = df_meth.loc[probe, mut_overlap].mean()
                wt_val = df_meth.loc[probe, wt_overlap].mean()
                report.append(f"    {probe}: Mutated={mut_val:.4f}, Wildtype={wt_val:.4f}, Delta={mut_val - wt_val:.4f}")
    else:
        report.append(f"  Insufficient overlap for methylation comparison (Mut={len(mut_overlap)}, WT={len(wt_overlap)})")
    
    # ---- LAYER 3: TRANSCRIPTOMICS (RNA-Seq) ----
    report.append(f"\n--- LAYER 3: TRANSCRIPTOMICS (RNA-Seq) ---")
    mut_rna = [p for p in mutated_patients if p in rna_patients]
    wt_rna = [p for p in wildtype_patients if p in rna_patients]
    
    if len(mut_rna) >= 2 and len(wt_rna) >= 2:
        report.append(f"  Comparing RNA: {len(mut_rna)} mutated vs {len(wt_rna)} wildtype patients")
        
        mean_mut_rna = df_rna[mut_rna].mean(axis=1)
        mean_wt_rna = df_rna[wt_rna].mean(axis=1)
        log2fc = np.log2((mean_mut_rna + 1) / (mean_wt_rna + 1))
        
        # Top upregulated genes in mutated group
        top_up = log2fc.sort_values(ascending=False).head(15)
        report.append(f"  Top 15 Upregulated Genes in {gene}-mutated:")
        for g, fc in top_up.items():
            report.append(f"    {g}: Log2FC = {fc:.2f}")
        
        # Top downregulated genes
        top_down = log2fc.sort_values(ascending=True).head(15)
        report.append(f"  Top 15 Downregulated Genes in {gene}-mutated:")
        for g, fc in top_down.items():
            report.append(f"    {g}: Log2FC = {fc:.2f}")
        
        # Check specific target genes
        target_rna_genes = ['MGMT', 'STAT3', 'EGFR', 'GRIA2', 'NLGN3', 'CD44', 'OLIG2', 'SOX2', 'CHI3L1', 'VIM']
        report.append(f"  Target Gene Expression in {gene}-mutated vs wildtype:")
        for tg in target_rna_genes:
            if tg in log2fc.index:
                fc = log2fc[tg]
                status = "UP" if fc > 1 else ("DOWN" if fc < -1 else "NO CHANGE")
                report.append(f"    {tg}: Log2FC = {fc:.2f} ({status})")
        
        # Compare to Normal Brain
        if gene in df_rna.index:
            tumor_expr = df_rna.loc[gene, mut_rna].mean() if len(mut_rna) > 0 else 0
            normal_expr = df_normal_cpm.loc[gene].mean() if gene in df_normal_cpm.index else 0
            report.append(f"  {gene} itself: Tumor Mean = {tumor_expr:.2f}, Normal Brain Mean = {normal_expr:.2f}")
    else:
        report.append(f"  Insufficient overlap for RNA comparison (Mut={len(mut_rna)}, WT={len(wt_rna)})")
    
    # ---- LAYER 4: PROTEOMICS ----
    report.append(f"\n--- LAYER 4: PROTEOMICS (Protein Abundance) ---")
    mut_prot = [p for p in mutated_patients if p in prot_patient_cols]
    wt_prot = [p for p in wildtype_patients if p in prot_patient_cols]
    
    if len(mut_prot) >= 2 and len(wt_prot) >= 2:
        report.append(f"  Comparing Proteomics: {len(mut_prot)} mutated vs {len(wt_prot)} wildtype patients")
        
        # Check specific target proteins
        target_proteins = ['STAT3', 'EGFR', 'PTEN', 'IDH1', 'GRIA2', 'OLIG2', 'CD44', 'MGMT', 'TP53']
        report.append(f"  Target Protein Abundance in {gene}-mutated vs wildtype:")
        for tp in target_proteins:
            tp_rows = df_prot[df_prot['Genes'] == tp]
            if not tp_rows.empty:
                row = tp_rows.iloc[0]
                mut_vals = pd.to_numeric(row[mut_prot], errors='coerce').dropna()
                wt_vals = pd.to_numeric(row[wt_prot], errors='coerce').dropna()
                if len(mut_vals) >= 2 and len(wt_vals) >= 2:
                    mut_mean = mut_vals.mean()
                    wt_mean = wt_vals.mean()
                    try:
                        stat, pval = mannwhitneyu(mut_vals, wt_vals, alternative='two-sided')
                        sig = "SIGNIFICANT" if pval < 0.05 else "not significant"
                        report.append(f"    {tp}: Mutated Mean={mut_mean:.2f}, WT Mean={wt_mean:.2f}, p={pval:.4f} ({sig})")
                    except:
                        report.append(f"    {tp}: Mutated Mean={mut_mean:.2f}, WT Mean={wt_mean:.2f}")
    else:
        report.append(f"  Insufficient overlap for proteomics comparison (Mut={len(mut_prot)}, WT={len(wt_prot)})")
    
    # ---- LAYER 5: PHOSPHOPROTEOMICS ----
    report.append(f"\n--- LAYER 5: PHOSPHOPROTEOMICS (Active Protein Switches) ---")
    mut_phos = [p for p in mutated_patients if p in phos_patient_cols]
    wt_phos = [p for p in wildtype_patients if p in phos_patient_cols]
    
    if len(mut_phos) >= 2 and len(wt_phos) >= 2:
        report.append(f"  Comparing Phosphoproteomics: {len(mut_phos)} mutated vs {len(wt_phos)} wildtype patients")
        
        target_phospho = ['STAT3_Y705', 'STAT3_S727', 'EGFR_S991', 'EGFR_T693', 'AKT1_S129', 'MAPK1_Y187', 'MAPK1_T185']
        report.append(f"  Target Phosphorylation Sites in {gene}-mutated vs wildtype:")
        for tp in target_phospho:
            tp_rows = df_phos[df_phos['PTM.CollapseKey'] == tp]
            if not tp_rows.empty:
                row = tp_rows.iloc[0]
                mut_vals = pd.to_numeric(row[mut_phos], errors='coerce').dropna()
                wt_vals = pd.to_numeric(row[wt_phos], errors='coerce').dropna()
                if len(mut_vals) >= 2 and len(wt_vals) >= 2:
                    mut_mean = mut_vals.mean()
                    wt_mean = wt_vals.mean()
                    try:
                        stat, pval = mannwhitneyu(mut_vals, wt_vals, alternative='two-sided')
                        sig = "SIGNIFICANT" if pval < 0.05 else "not significant"
                        report.append(f"    {tp}: Mutated Mean={mut_mean:.2f}, WT Mean={wt_mean:.2f}, p={pval:.4f} ({sig})")
                    except:
                        report.append(f"    {tp}: Mutated Mean={mut_mean:.2f}, WT Mean={wt_mean:.2f}")
    else:
        report.append(f"  Insufficient overlap for phosphoproteomics comparison (Mut={len(mut_phos)}, WT={len(wt_phos)})")

# =========================================================
# TREATMENT EFFECT ANALYSIS (TMZ Treated vs Untreated)
# =========================================================
report.append("\n\n" + "=" * 80)
report.append("TREATMENT EFFECT ANALYSIS")
report.append("Comparing TMZ-Treated vs Untreated across all omic layers")
report.append("=" * 80)

chemo_col = 'Chemo_status (TMZ treated=1;un-treated=0)'
# Get treated vs untreated patient IDs from WES clinical
treated_ids = df_wes_clin[df_wes_clin[chemo_col] == 1].index.tolist()
untreated_ids = df_wes_clin[df_wes_clin[chemo_col] == 0].index.tolist()

report.append(f"TMZ Treated: {len(treated_ids)} patients")
report.append(f"TMZ Untreated: {len(untreated_ids)} patients")

# Treatment effect on RNA
treated_rna = [p for p in treated_ids if p in rna_patients]
untreated_rna = [p for p in untreated_ids if p in rna_patients]

if len(treated_rna) >= 2 and len(untreated_rna) >= 2:
    report.append(f"\n--- RNA Expression: TMZ Treated ({len(treated_rna)}) vs Untreated ({len(untreated_rna)}) ---")
    mean_treated = df_rna[treated_rna].mean(axis=1)
    mean_untreated = df_rna[untreated_rna].mean(axis=1)
    log2fc_treat = np.log2((mean_treated + 1) / (mean_untreated + 1))
    
    target_rna_genes = ['MGMT', 'STAT3', 'EGFR', 'GRIA2', 'NLGN3', 'CD44', 'OLIG2', 'SOX2', 'CHI3L1', 'VIM', 'PTEN', 'TP53', 'ATRX']
    for tg in target_rna_genes:
        if tg in log2fc_treat.index:
            fc = log2fc_treat[tg]
            status = "UP in Treated" if fc > 0.5 else ("DOWN in Treated" if fc < -0.5 else "NO CHANGE")
            report.append(f"  {tg}: Log2FC = {fc:.2f} ({status})")

# Treatment effect on Methylation
treated_meth = [p for p in treated_ids if p in meth_patients]
untreated_meth = [p for p in untreated_ids if p in meth_patients]

if len(treated_meth) >= 2 and len(untreated_meth) >= 2:
    report.append(f"\n--- Methylation: TMZ Treated ({len(treated_meth)}) vs Untreated ({len(untreated_meth)}) ---")
    mean_t_meth = df_meth[treated_meth].mean(axis=1)
    mean_u_meth = df_meth[untreated_meth].mean(axis=1)
    delta_treat_meth = mean_t_meth - mean_u_meth
    
    report.append("  Top 10 probes with higher methylation in Treated:")
    for probe, delta in delta_treat_meth.sort_values(ascending=False).head(10).items():
        report.append(f"    {probe}: Delta Beta = {delta:.4f}")
    
    report.append("  Top 10 probes with lower methylation in Treated:")
    for probe, delta in delta_treat_meth.sort_values(ascending=True).head(10).items():
        report.append(f"    {probe}: Delta Beta = {delta:.4f}")

# Treatment effect on Proteomics
treated_prot = [p for p in treated_ids if p in prot_patient_cols]
untreated_prot = [p for p in untreated_ids if p in prot_patient_cols]

if len(treated_prot) >= 2 and len(untreated_prot) >= 2:
    report.append(f"\n--- Proteomics: TMZ Treated ({len(treated_prot)}) vs Untreated ({len(untreated_prot)}) ---")
    target_proteins = ['STAT3', 'EGFR', 'PTEN', 'IDH1', 'GRIA2', 'OLIG2', 'CD44', 'MGMT', 'TP53']
    for tp in target_proteins:
        tp_rows = df_prot[df_prot['Genes'] == tp]
        if not tp_rows.empty:
            row = tp_rows.iloc[0]
            t_vals = pd.to_numeric(row[treated_prot], errors='coerce').dropna()
            u_vals = pd.to_numeric(row[untreated_prot], errors='coerce').dropna()
            if len(t_vals) >= 2 and len(u_vals) >= 2:
                t_mean = t_vals.mean()
                u_mean = u_vals.mean()
                try:
                    stat, pval = mannwhitneyu(t_vals, u_vals, alternative='two-sided')
                    sig = "SIGNIFICANT" if pval < 0.05 else "not significant"
                    report.append(f"  {tp}: Treated Mean={t_mean:.2f}, Untreated Mean={u_mean:.2f}, p={pval:.4f} ({sig})")
                except:
                    report.append(f"  {tp}: Treated Mean={t_mean:.2f}, Untreated Mean={u_mean:.2f}")

# =========================================================
# CO-MUTATION PATTERNS
# =========================================================
report.append("\n\n" + "=" * 80)
report.append("CO-MUTATION PATTERNS")
report.append("Which mutations tend to occur together?")
report.append("=" * 80)

for g1 in ['IDH1', 'TP53', 'ATRX', 'PTEN', 'EGFR', 'CIC']:
    for g2 in ['IDH1', 'TP53', 'ATRX', 'PTEN', 'EGFR', 'CIC']:
        if g1 >= g2:
            continue
        if g1 in df_mut_bin.columns and g2 in df_mut_bin.columns:
            both = ((df_mut_bin[g1] == 1) & (df_mut_bin[g2] == 1)).sum()
            g1_only = ((df_mut_bin[g1] == 1) & (df_mut_bin[g2] == 0)).sum()
            g2_only = ((df_mut_bin[g1] == 0) & (df_mut_bin[g2] == 1)).sum()
            neither = ((df_mut_bin[g1] == 0) & (df_mut_bin[g2] == 0)).sum()
            total = len(df_mut_bin)
            co_rate = both / total * 100
            if co_rate > 2:
                report.append(f"  {g1} + {g2}: Co-mutated in {both} patients ({co_rate:.1f}%)")

# Save report
report_text = "\n".join(report)
with open(output_file, 'w') as f:
    f.write(report_text)

print(f"\nReport written to: {output_file}")
print(f"Report length: {len(report)} lines")
