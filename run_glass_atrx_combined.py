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

# 1. Load RNA-seq
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

log(f"Matched RNA patients: {len(matched_tp_r1)}")

# 2. Load mutations using chunking to handle truncation
log("\nLoading coding mutations from truncated variants.anno.csv.gz...")
df_coding = pd.DataFrame()
try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.anno.csv.gz"), chunksize=200000)
    for chunk in reader:
        coding_chunk = chunk[chunk['variant_classification'].isin([
            'MISSENSE', 'NONSENSE', 'FRAME_SHIFT_DEL', 'FRAME_SHIFT_INS',
            'SPLICE_SITE', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'NONSTOP'
        ])][['variant_id', 'gene_symbol', 'variant_classification']]
        df_coding = pd.concat([df_coding, coding_chunk], ignore_index=True)
except EOFError:
    pass

log(f"Recovered coding mutations: {len(df_coding)}")

log("\nLoading genotypes from truncated variants.passgeno.csv.gz...")
df_geno = pd.DataFrame()
try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=200000)
    for chunk in reader:
        pass_chunk = chunk[chunk['ssm2_pass_call'] == 1][['case_barcode', 'variant_id']]
        df_geno = pd.concat([df_geno, pass_chunk], ignore_index=True)
except EOFError:
    pass

log(f"Recovered genotypes: {len(df_geno)}")

# Merge in memory
if not df_coding.empty and not df_geno.empty:
    df_mut = pd.merge(df_geno, df_coding, on='variant_id')
    log(f"Merged mutations: {len(df_mut)}")
    
    # Get ATRX mutants
    atrx_mut_cases = set(df_mut[df_mut['gene_symbol'] == 'ATRX']['case_barcode'].unique())
    log(f"ATRX mutant cases in genotype data: {len(atrx_mut_cases)}")

    # Map patients
    atrx_mut_matched = []
    atrx_wt_matched = []
    for p in matched_tp_r1:
        var_p = p.replace(".", "-")
        if var_p in atrx_mut_cases:
            atrx_mut_matched.append(p)
        else:
            # Check if patient exists in our genotype dataset to make sure they were sequenced
            if var_p in df_mut['case_barcode'].unique():
                atrx_wt_matched.append(p)

    log(f"Matched TP->R1 with ATRX-Mut: {len(atrx_mut_matched)}")
    log(f"Matched TP->R1 with ATRX-WT: {len(atrx_wt_matched)}")

    key_markers = ['CHI3L1','CD44','STAT3','GRIA2','NLGN3','EGFR','VIM','CD163','AIF1','OLIG2','SOX2','MKI67']

    log(f"\n{'Gene':<12} {'ATRX-Mut TP':>12} {'ATRX-Mut R1':>12} {'Mut FC':>8} {'ATRX-WT TP':>12} {'ATRX-WT R1':>12} {'WT FC':>8}")
    log("-"*85)

    for g in key_markers:
        if g in df_rna.index:
            mut_tp = df_rna.loc[g, [tp_samples[p] for p in atrx_mut_matched]].median() if len(atrx_mut_matched) > 0 else 0
            mut_r1 = df_rna.loc[g, [r1_samples[p] for p in atrx_mut_matched]].median() if len(atrx_mut_matched) > 0 else 0
            mut_fc = np.log2((mut_r1+1)/(mut_tp+1)) if len(atrx_mut_matched) > 0 else 0
            
            wt_tp = df_rna.loc[g, [tp_samples[p] for p in atrx_wt_matched]].median() if len(atrx_wt_matched) > 0 else 0
            wt_r1 = df_rna.loc[g, [r1_samples[p] for p in atrx_wt_matched]].median() if len(atrx_wt_matched) > 0 else 0
            wt_fc = np.log2((wt_r1+1)/(wt_tp+1)) if len(atrx_wt_matched) > 0 else 0
            
            log(f"{g:<12} {mut_tp:>12.2f} {mut_r1:>12.2f} {mut_fc:>8.2f} {wt_tp:>12.2f} {wt_r1:>12.2f} {wt_fc:>8.2f}")
else:
    log("Failed to load mutation/genotype data.")

with open(out, 'w') as f:
    f.write("\n".join(report))
log(f"\nResults saved to {out}")
