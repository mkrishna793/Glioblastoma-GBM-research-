import pandas as pd
import numpy as np
import gzip
import os, sys
sys.stdout.reconfigure(encoding='utf-8')

data_dir = r"D:\new data of the GBM research"
out = r"D:\new data of the GBM research\GLASS_longitudinal_report.txt"
report = []

def log(msg):
    report.append(msg)
    print(msg)

log("="*80)
log("GLASS LONGITUDINAL TUMOR EVOLUTION ANALYSIS")
log("="*80)

# =====================================================================
# STEP 1: LOAD RNA-SEQ AND BUILD MATCHED PATIENT PAIRS
# =====================================================================
log("\n[STEP 1] Loading RNA-seq TPM matrix...")
df_rna = pd.read_csv(os.path.join(data_dir, "gene_tpm_matrix_all_samples.tsv"), sep="\t", index_col=0)
df_rna.columns = df_rna.columns.str.strip()

# Parse barcodes into patient ID + timepoint
sample_meta = []
for col in df_rna.columns:
    parts = col.replace("-", ".").split(".")
    pid = ".".join(parts[:3])
    tp = parts[3] if len(parts) > 3 else "UNK"
    sample_meta.append({"sample": col, "patient": pid, "timepoint": tp})

meta = pd.DataFrame(sample_meta)
tp_samples = meta[meta["timepoint"] == "TP"].drop_duplicates(subset="patient").set_index("patient")["sample"]
r1_samples = meta[meta["timepoint"] == "R1"].drop_duplicates(subset="patient").set_index("patient")["sample"]
r2_samples = meta[meta["timepoint"] == "R2"].drop_duplicates(subset="patient").set_index("patient")["sample"]

matched_tp_r1 = sorted(set(tp_samples.index) & set(r1_samples.index))
matched_tp_r1_r2 = sorted(set(matched_tp_r1) & set(r2_samples.index))

log(f"  Total RNA samples: {len(df_rna.columns)}")
log(f"  Primary (TP): {len(tp_samples)}, Recurrence (R1): {len(r1_samples)}, 2nd Recurrence (R2): {len(r2_samples)}")
log(f"  Matched TP→R1 patients: {len(matched_tp_r1)}")
log(f"  Matched TP→R1→R2 patients: {len(matched_tp_r1_r2)}")

# =====================================================================
# STEP 2: TRACE KEY GENE EXPRESSION CHANGES (TP → R1)
# =====================================================================
log("\n" + "="*80)
log("[STEP 2] GENE EXPRESSION EVOLUTION: PRIMARY → RECURRENCE (N=%d matched)" % len(matched_tp_r1))
log("="*80)

targets = ['CHI3L1','CD44','STAT3','EGFR','GRIA2','NLGN3','OLIG2','SOX2',
           'VIM','MGMT','MKI67','PDGFRA','NF1','PTEN','CD163','AIF1',
           'TGFB1','IL6','TNF','HIF1A','VEGFA','MMP2','MMP9','FN1',
           'CXCL12','CCL2','CSF1','PTPRC','CD3E','CD8A','FOXP3','CD274']

tp_expr = df_rna[[tp_samples[p] for p in matched_tp_r1]]
tp_expr.columns = matched_tp_r1
r1_expr = df_rna[[r1_samples[p] for p in matched_tp_r1]]
r1_expr.columns = matched_tp_r1

log("\n--- Per-Gene Median Change (Primary → Recurrence) ---")
log(f"{'Gene':<12} {'TP Median':>12} {'R1 Median':>12} {'Log2FC':>10} {'Direction':<20}")
log("-"*70)

for g in targets:
    if g in tp_expr.index:
        tp_vals = tp_expr.loc[g]
        r1_vals = r1_expr.loc[g]
        tp_med = tp_vals.median()
        r1_med = r1_vals.median()
        fc = np.log2((r1_med + 1) / (tp_med + 1))
        direction = "UP in recurrence" if fc > 0.5 else ("DOWN in recurrence" if fc < -0.5 else "STABLE")
        log(f"{g:<12} {tp_med:>12.2f} {r1_med:>12.2f} {fc:>10.2f} {direction:<20}")

# =====================================================================
# STEP 3: TRACE TP → R1 → R2 EVOLUTION (TRIPLE-MATCHED)
# =====================================================================
log("\n" + "="*80)
log("[STEP 3] TRIPLE-TIMEPOINT EVOLUTION: TP → R1 → R2 (N=%d patients)" % len(matched_tp_r1_r2))
log("="*80)

if len(matched_tp_r1_r2) > 0:
    log(f"\n{'Gene':<12} {'TP Median':>12} {'R1 Median':>12} {'R2 Median':>12} {'FC(TP→R1)':>10} {'FC(R1→R2)':>10} {'FC(TP→R2)':>10}")
    log("-"*85)
    
    tp2 = df_rna[[tp_samples[p] for p in matched_tp_r1_r2]]
    r12 = df_rna[[r1_samples[p] for p in matched_tp_r1_r2]]
    r22 = df_rna[[r2_samples[p] for p in matched_tp_r1_r2]]
    
    for g in targets:
        if g in tp2.index:
            t = tp2.loc[g].median()
            r1 = r12.loc[g].median()
            r2 = r22.loc[g].median()
            fc1 = np.log2((r1+1)/(t+1))
            fc2 = np.log2((r2+1)/(r1+1))
            fc3 = np.log2((r2+1)/(t+1))
            log(f"{g:<12} {t:>12.2f} {r1:>12.2f} {r2:>12.2f} {fc1:>10.2f} {fc2:>10.2f} {fc3:>10.2f}")
else:
    log("  No triple-matched patients found.")

# =====================================================================
# STEP 4: CD44/YKL-40 DECONVOLUTION VALIDATION
# =====================================================================
log("\n" + "="*80)
log("[STEP 4] CD44 vs YKL-40 DECONVOLUTION VALIDATION")
log("="*80)

# Use macrophage markers (CD163, AIF1/Iba1) as proxy for immune infiltration
# If CD44 correlates with macrophage markers but CHI3L1 does not, our finding is validated
macro_markers = ['CD163', 'AIF1', 'CSF1R', 'ITGAM']
cancer_markers = ['CHI3L1', 'SOX2', 'OLIG2', 'NESTIN']

log("\n--- Correlation of CD44 and CHI3L1 with Macrophage Markers (All samples) ---")
all_samples = df_rna.columns.tolist()

for marker in macro_markers:
    if marker in df_rna.index and 'CD44' in df_rna.index:
        cd44_corr = df_rna.loc['CD44', all_samples].corr(df_rna.loc[marker, all_samples])
        chi3l1_corr = df_rna.loc['CHI3L1', all_samples].corr(df_rna.loc[marker, all_samples]) if 'CHI3L1' in df_rna.index else float('nan')
        log(f"  {marker} vs CD44:   r = {cd44_corr:.4f}")
        log(f"  {marker} vs CHI3L1: r = {chi3l1_corr:.4f}")
        log("")

# Also check CD44 vs CHI3L1 direct correlation
if 'CD44' in df_rna.index and 'CHI3L1' in df_rna.index:
    direct = df_rna.loc['CD44', all_samples].corr(df_rna.loc['CHI3L1', all_samples])
    log(f"  Direct CD44 vs CHI3L1 correlation: r = {direct:.4f}")

# Per-patient change analysis: does CD44 change track with macrophage change?
log("\n--- Per-Patient Recurrence Change Correlation ---")
if len(matched_tp_r1) > 5:
    delta_cd44 = []
    delta_chi3l1 = []
    delta_cd163 = []
    delta_aif1 = []
    
    for p in matched_tp_r1:
        tp_s = tp_samples[p]
        r1_s = r1_samples[p]
        if 'CD44' in df_rna.index:
            delta_cd44.append(df_rna.loc['CD44', r1_s] - df_rna.loc['CD44', tp_s])
        if 'CHI3L1' in df_rna.index:
            delta_chi3l1.append(df_rna.loc['CHI3L1', r1_s] - df_rna.loc['CHI3L1', tp_s])
        if 'CD163' in df_rna.index:
            delta_cd163.append(df_rna.loc['CD163', r1_s] - df_rna.loc['CD163', tp_s])
        if 'AIF1' in df_rna.index:
            delta_aif1.append(df_rna.loc['AIF1', r1_s] - df_rna.loc['AIF1', tp_s])
    
    delta_cd44 = pd.Series(delta_cd44)
    delta_chi3l1 = pd.Series(delta_chi3l1)
    delta_cd163 = pd.Series(delta_cd163)
    delta_aif1 = pd.Series(delta_aif1)
    
    log(f"  Delta(CD44) vs Delta(CD163):   r = {delta_cd44.corr(delta_cd163):.4f}")
    log(f"  Delta(CD44) vs Delta(AIF1):    r = {delta_cd44.corr(delta_aif1):.4f}")
    log(f"  Delta(CHI3L1) vs Delta(CD163): r = {delta_chi3l1.corr(delta_cd163):.4f}")
    log(f"  Delta(CHI3L1) vs Delta(AIF1):  r = {delta_chi3l1.corr(delta_aif1):.4f}")
    log(f"  Delta(CD44) vs Delta(CHI3L1):  r = {delta_cd44.corr(delta_chi3l1):.4f}")

# =====================================================================
# STEP 5: LOAD MUTATIONS - IDENTIFY IDH/ATRX/TP53 STATUS
# =====================================================================
log("\n" + "="*80)
log("[STEP 5] MUTATION STATUS FROM GLASS VARIANTS")
log("="*80)

log("  Loading annotated variants...")
with gzip.open(os.path.join(data_dir, "variants.anno.csv.gz"), "rt") as f:
    df_var_anno = pd.read_csv(f)

# Filter for coding mutations
coding = df_var_anno[df_var_anno['variant_classification'].isin([
    'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
    'Splice_Site', 'In_Frame_Del', 'In_Frame_Ins', 'Nonstop_Mutation'
])]

log(f"  Total annotated variants: {len(df_var_anno)}")
log(f"  Coding variants: {len(coding)}")

# Get variant IDs for key genes
key_genes = ['IDH1', 'TP53', 'ATRX', 'PTEN', 'EGFR', 'NF1', 'CIC', 'PIK3CA', 'PDGFRA']
for g in key_genes:
    g_vars = coding[coding['gene_symbol'] == g]
    log(f"  {g}: {len(g_vars)} coding variants found")

# Load pass genotypes to map variants to patients
log("\n  Loading pass genotypes...")
with gzip.open(os.path.join(data_dir, "variants.passgeno.csv.gz"), "rt") as f:
    df_geno = pd.read_csv(f)

log(f"  Total genotype calls: {len(df_geno)}")

# Get IDH1 variant IDs
idh1_var_ids = coding[coding['gene_symbol'] == 'IDH1']['variant_id'].tolist()
atrx_var_ids = coding[coding['gene_symbol'] == 'ATRX']['variant_id'].tolist()
tp53_var_ids = coding[coding['gene_symbol'] == 'TP53']['variant_id'].tolist()

# Find patients with these mutations (use pass calls)
idh1_patients = set(df_geno[df_geno['variant_id'].isin(idh1_var_ids) & (df_geno['ssm2_pass_call'] == 1)]['case_barcode'].unique())
atrx_patients = set(df_geno[df_geno['variant_id'].isin(atrx_var_ids) & (df_geno['ssm2_pass_call'] == 1)]['case_barcode'].unique())
tp53_patients = set(df_geno[df_geno['variant_id'].isin(tp53_var_ids) & (df_geno['ssm2_pass_call'] == 1)]['case_barcode'].unique())

log(f"\n  IDH1-mutant patients: {len(idh1_patients)}")
log(f"  ATRX-mutant patients: {len(atrx_patients)}")
log(f"  TP53-mutant patients: {len(tp53_patients)}")

# =====================================================================
# STEP 6: ATRX-STRATIFIED RECURRENCE ANALYSIS
# =====================================================================
log("\n" + "="*80)
log("[STEP 6] ATRX-STRATIFIED TUMOR EVOLUTION AT RECURRENCE")
log("="*80)

# Map GLASS patient barcodes to RNA patient IDs
# RNA patients use "." separator, variant patients use "-"
rna_to_var = {}
for p in matched_tp_r1:
    var_p = p.replace(".", "-")
    rna_to_var[p] = var_p

atrx_mut_matched = [p for p in matched_tp_r1 if rna_to_var[p] in atrx_patients]
atrx_wt_matched = [p for p in matched_tp_r1 if rna_to_var[p] not in atrx_patients]

log(f"  Matched TP→R1 with ATRX-Mut: {len(atrx_mut_matched)}")
log(f"  Matched TP→R1 with ATRX-WT: {len(atrx_wt_matched)}")

key_markers = ['CHI3L1','CD44','STAT3','GRIA2','NLGN3','EGFR','VIM','CD163','AIF1','OLIG2','SOX2','MKI67']

if len(atrx_mut_matched) >= 3 and len(atrx_wt_matched) >= 3:
    log(f"\n{'Gene':<12} {'ATRX-Mut TP':>12} {'ATRX-Mut R1':>12} {'Mut FC':>8} {'ATRX-WT TP':>12} {'ATRX-WT R1':>12} {'WT FC':>8}")
    log("-"*80)
    
    for g in key_markers:
        if g in df_rna.index:
            # ATRX-Mut
            mut_tp = df_rna.loc[g, [tp_samples[p] for p in atrx_mut_matched]].median()
            mut_r1 = df_rna.loc[g, [r1_samples[p] for p in atrx_mut_matched]].median()
            mut_fc = np.log2((mut_r1+1)/(mut_tp+1))
            # ATRX-WT
            wt_tp = df_rna.loc[g, [tp_samples[p] for p in atrx_wt_matched]].median()
            wt_r1 = df_rna.loc[g, [r1_samples[p] for p in atrx_wt_matched]].median()
            wt_fc = np.log2((wt_r1+1)/(wt_tp+1))
            log(f"{g:<12} {mut_tp:>12.2f} {mut_r1:>12.2f} {mut_fc:>8.2f} {wt_tp:>12.2f} {wt_r1:>12.2f} {wt_fc:>8.2f}")

# =====================================================================
# STEP 7: METHYLATION EVOLUTION (TP → R1)
# =====================================================================
log("\n" + "="*80)
log("[STEP 7] METHYLATION EVOLUTION: PRIMARY → RECURRENCE")
log("="*80)

log("  Loading 450K methylation...")
df_meth = pd.read_csv(os.path.join(data_dir, "beta.450.tsv"), sep="\t", index_col=0)
df_meth.columns = df_meth.columns.str.strip()

# Parse methylation barcodes
meth_meta = []
for col in df_meth.columns:
    parts = col.replace("-", ".").split(".")
    # TCGA.DH.A669.TP.12M.450.QDG6WF
    pid = ".".join(parts[:3])
    tp = parts[3] if len(parts) > 3 else "UNK"
    meth_meta.append({"sample": col, "patient": pid, "timepoint": tp})

meth_meta = pd.DataFrame(meth_meta)
meth_tp = meth_meta[meth_meta["timepoint"] == "TP"].set_index("patient")["sample"]
meth_r1 = meth_meta[meth_meta["timepoint"] == "R1"].set_index("patient")["sample"]
meth_matched = sorted(set(meth_tp.index) & set(meth_r1.index))

log(f"  Matched TP→R1 methylation patients: {len(meth_matched)}")

# MGMT probes
mgmt_probes = ['cg12434587', 'cg12981137']
log("\n--- MGMT Promoter Methylation Evolution ---")
for probe in mgmt_probes:
    if probe in df_meth.index and len(meth_matched) > 0:
        tp_vals = df_meth.loc[probe, [meth_tp[p] for p in meth_matched if p in meth_tp.index]]
        r1_vals = df_meth.loc[probe, [meth_r1[p] for p in meth_matched if p in meth_r1.index]]
        log(f"  {probe}: TP Median = {tp_vals.median():.4f}, R1 Median = {r1_vals.median():.4f}, Delta = {r1_vals.median() - tp_vals.median():.4f}")

# Global methylation shift
log("\n--- Global Methylation Shift (Top changes) ---")
if len(meth_matched) > 0:
    tp_mean = df_meth[[meth_tp[p] for p in meth_matched if p in meth_tp.index]].mean(axis=1)
    r1_mean = df_meth[[meth_r1[p] for p in meth_matched if p in meth_r1.index]].mean(axis=1)
    delta = r1_mean - tp_mean
    
    log("  Top 10 Hypermethylated at Recurrence (genes locked):")
    for probe, d in delta.sort_values(ascending=False).head(10).items():
        log(f"    {probe}: Delta = {d:.4f}")
    
    log("  Top 10 Hypomethylated at Recurrence (genes unlocked):")
    for probe, d in delta.sort_values(ascending=True).head(10).items():
        log(f"    {probe}: Delta = {d:.4f}")

# =====================================================================
# STEP 8: COPY NUMBER CHANGES AT RECURRENCE
# =====================================================================
log("\n" + "="*80)
log("[STEP 8] COPY NUMBER ALTERATIONS")
log("="*80)

log("  Loading gene copy number data...")
with gzip.open(os.path.join(data_dir, "variants.gene_copy_number.csv.gz"), "rt") as f:
    df_cn = pd.read_csv(f)

log(f"  Total CN calls: {len(df_cn)}")

# Check key gene amplifications/deletions
for g in ['EGFR', 'PDGFRA', 'CDK4', 'MDM2', 'PTEN', 'CDKN2A', 'RB1', 'NF1']:
    g_cn = df_cn[df_cn['gene_symbol'] == g]
    amp = len(g_cn[g_cn['hlvl_call'] >= 2])
    dele = len(g_cn[g_cn['hlvl_call'] <= -2])
    log(f"  {g}: Amplified in {amp} samples, Deleted in {dele} samples")

# =====================================================================
# SAVE REPORT
# =====================================================================
with open(out, 'w') as f:
    f.write("\n".join(report))

log("\n\nReport saved to: " + out)
