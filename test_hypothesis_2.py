import pandas as pd
import numpy as np
import os
from scipy.stats import wilcoxon

data_dir = r"D:\new data of the GBM research"
out_path = r"D:\research of the GBM\hypothesis_2_results.txt"

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

# Define subtype markers
signatures = {
    'NEU': ['GRIA2', 'NLGN3', 'SYT1', 'SNAP25'],
    'AFM': ['STAT3', 'CHI3L1', 'VIM', 'GFAP'],
    'PPR': ['OLIG2', 'SOX2', 'PROM1', 'NES'],
    'IME': ['CD44', 'CD163', 'CSF1R', 'ITGAM']
}

# Verify all markers exist in RNA-seq matrix
for sub, markers in signatures.items():
    existing = [m for m in markers if m in df_rna.index]
    print(f"  {sub} markers found: {len(existing)} / {len(markers)} ({existing})")
    signatures[sub] = existing

# Calculate log2(TPM + 1)
df_log = np.log2(df_rna + 1)

# Calculate raw signature scores (median log2 expression of marker genes)
raw_scores = {}
for sub, markers in signatures.items():
    raw_scores[sub] = df_log.loc[markers].median(axis=0)

df_scores = pd.DataFrame(raw_scores)

# Standardize scores across all samples using z-score (to make them comparable)
df_zscores = (df_scores - df_scores.mean()) / df_scores.std()

# Map samples to patient and timepoint
df_zscores = df_zscores.join(df_rna_meta.set_index('sample'))

# Separate TP and R1 scores for matched patients
scores_tp = df_zscores[(df_zscores['patient_id'].isin(matched_rna_patients)) & (df_zscores['timepoint'] == 'TP')]
scores_r1 = df_zscores[(df_zscores['patient_id'].isin(matched_rna_patients)) & (df_zscores['timepoint'] == 'R1')]

# Group by patient_id to resolve duplicate samples per patient (take the mean)
scores_tp = scores_tp.groupby('patient_id')[['NEU', 'AFM', 'PPR', 'IME']].mean()
scores_r1 = scores_r1.groupby('patient_id')[['NEU', 'AFM', 'PPR', 'IME']].mean()

# Reindex both to match patient order
scores_tp = scores_tp.reindex(matched_rna_patients)
scores_r1 = scores_r1.reindex(matched_rna_patients)

# Perform Wilcoxon signed-rank tests for each subtype
with open(out_path, "w") as f:
    f.write("="*80 + "\n")
    f.write("HYPOTHESIS 2 TESTING: LONGITUDINAL PHENOTYPIC SUBTYPE DRIFT (TP -> R1)\n")
    f.write("="*80 + "\n\n")
    
    f.write(f"Matched Patients analyzed: {len(matched_rna_patients)}\n\n")
    f.write(f"{'Subtype':<10} | {'TP Mean Z':>10} | {'R1 Mean Z':>10} | {'Wilcoxon stat':>15} | {'p-value':>12} | {'Direction':<20}\n")
    f.write("-"*85 + "\n")
    
    for sub in ['NEU', 'AFM', 'PPR', 'IME']:
        tp_vals = scores_tp[sub].values
        r1_vals = scores_r1[sub].values
        
        tp_mean = np.mean(tp_vals)
        r1_mean = np.mean(r1_vals)
        
        stat, p_val = wilcoxon(r1_vals, tp_vals)
        
        direction = "INCREASES in R1" if r1_mean > tp_mean else "DECREASES in R1"
        direction_f = f"{direction} (significant)" if p_val < 0.05 else f"{direction} (non-sig)"
        
        f.write(f"{sub:<10} | {tp_mean:>10.4f} | {r1_mean:>10.4f} | {stat:>15.1f} | {p_val:>12.4e} | {direction_f:<20}\n")
        
    # Classify patient subtype switching
    # A patient is assigned the subtype with the highest z-score in that sample
    tp_assigned = scores_tp[['NEU', 'AFM', 'PPR', 'IME']].idxmax(axis=1)
    r1_assigned = scores_r1[['NEU', 'AFM', 'PPR', 'IME']].idxmax(axis=1)
    
    df_switch = pd.DataFrame({'TP': tp_assigned, 'R1': r1_assigned})
    f.write("\n\nSubtype Transition Matrix (Number of Patients):\n")
    f.write(pd.crosstab(df_switch['TP'], df_switch['R1']).to_string() + "\n\n")
    
    # Calculate overall switch rate
    switched = (df_switch['TP'] != df_switch['R1']).sum()
    f.write(f"Total patients analyzed: {len(df_switch)}\n")
    f.write(f"Patients who switched subtype: {switched} ({switched/len(df_switch)*100:.2f}%)\n")
    f.write(f"Patients who remained stable: {len(df_switch)-switched} ({(len(df_switch)-switched)/len(df_switch)*100:.2f}%)\n")

print(f"Results saved to {out_path}")
# Print summary to stdout
print("\nSubtype Transition Summary:")
df_switch = pd.DataFrame({'TP': scores_tp[['NEU', 'AFM', 'PPR', 'IME']].idxmax(axis=1),
                          'R1': scores_r1[['NEU', 'AFM', 'PPR', 'IME']].idxmax(axis=1)})
print(pd.crosstab(df_switch['TP'], df_switch['R1']))
