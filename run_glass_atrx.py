import pandas as pd
import numpy as np
import os, sys

sys.stdout.reconfigure(encoding='utf-8')
data_dir = r"D:\new data of the GBM research"
out = r"D:\new data of the GBM research\GLASS_ATRX_stratified_results.txt"
report = []

def log(msg):
    report.append(msg)
    print(msg)

log("="*80)
log("ATRX-STRATIFIED LONGITUDINAL RNA EVOLUTION (GLASS)")
log("="*80)

# Load RNA-seq
log("Loading RNA-seq...")
df_rna = pd.read_csv(os.path.join(data_dir, "gene_tpm_matrix_all_samples.tsv"), sep="\t", index_col=0)
df_rna.columns = df_rna.columns.str.strip()

# Parse barcodes
sample_meta = []
for col in df_rna.columns:
    parts = col.replace("-", ".").split(".")
    pid = ".".join(parts[:3])
    tp = parts[3] if len(parts) > 3 else "UNK"
    sample_meta.append({"sample": col, "patient": pid, "timepoint": tp})

meta = pd.DataFrame(sample_meta)
tp_samples = meta[meta["timepoint"] == "TP"].drop_duplicates(subset="patient").set_index("patient")["sample"]
r1_samples = meta[meta["timepoint"] == "R1"].drop_duplicates(subset="patient").set_index("patient")["sample"]
matched_tp_r1 = sorted(set(tp_samples.index) & set(r1_samples.index))

log(f"Matched patients: {len(matched_tp_r1)}")

# Load salvaged mutations
log("Loading salvaged mutations...")
df_mut = pd.read_csv(os.path.join(data_dir, "GLASS_recovered_mutations.txt"), sep="\t")

# Get ATRX mutants
atrx_mut_cases = set(df_mut[df_mut['gene_symbol'] == 'ATRX']['case_barcode'].unique())
print(f"ATRX mutant cases in database: {len(atrx_mut_cases)}")

# Align cases (convert rna case to variant case: replace . with -)
atrx_mut_matched = []
atrx_wt_matched = []
for p in matched_tp_r1:
    var_p = p.replace(".", "-")
    if var_p in atrx_mut_cases:
        atrx_mut_matched.append(p)
    else:
        # Check if the case exists in df_mut to make sure it was sequenced
        if var_p in df_mut['case_barcode'].unique():
            atrx_wt_matched.append(p)

log(f"Matched TP->R1 with ATRX-Mut: {len(atrx_mut_matched)}")
log(f"Matched TP->R1 with ATRX-WT: {len(atrx_wt_matched)}")

key_markers = ['CHI3L1','CD44','STAT3','GRIA2','NLGN3','EGFR','VIM','CD163','AIF1','OLIG2','SOX2','MKI67']

log(f"\n{'Gene':<12} {'ATRX-Mut TP':>12} {'ATRX-Mut R1':>12} {'Mut FC':>8} {'ATRX-WT TP':>12} {'ATRX-WT R1':>12} {'WT FC':>8}")
log("-"*85)

for g in key_markers:
    if g in df_rna.index:
        # ATRX-Mut
        mut_tp_vals = df_rna.loc[g, [tp_samples[p] for p in atrx_mut_matched]].median() if len(atrx_mut_matched) > 0 else 0
        mut_r1_vals = df_rna.loc[g, [r1_samples[p] for p in atrx_mut_matched]].median() if len(atrx_mut_matched) > 0 else 0
        mut_fc = np.log2((mut_r1_vals+1)/(mut_tp_vals+1)) if len(atrx_mut_matched) > 0 else 0
        # ATRX-WT
        wt_tp_vals = df_rna.loc[g, [tp_samples[p] for p in atrx_wt_matched]].median() if len(atrx_wt_matched) > 0 else 0
        wt_r1_vals = df_rna.loc[g, [r1_samples[p] for p in atrx_wt_matched]].median() if len(atrx_wt_matched) > 0 else 0
        wt_fc = np.log2((wt_r1_vals+1)/(wt_tp_vals+1)) if len(atrx_wt_matched) > 0 else 0
        
        log(f"{g:<12} {mut_tp_vals:>12.2f} {mut_r1_vals:>12.2f} {mut_fc:>8.2f} {wt_tp_vals:>12.2f} {wt_r1_vals:>12.2f} {wt_fc:>8.2f}")

with open(out, 'w') as f:
    f.write("\n".join(report))
