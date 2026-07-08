"""Stage 5: Epigenome-to-Expression Check.
Tests promoter methylation change (450K) vs matched RNA expression change (TP->R1).
"""
import pandas as pd, numpy as np, os, gc, sys
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8')

data_dir = r"D:\new data of the GBM research"
out_dir  = r"D:\research of the GBM"

# ── STEP 1: Load matched probes dynamically from cpg_matched_probes.txt ──
print("STEP 1: Loading matched probes from local reference...")
chi3l1_probes = []
cd44_probes = []
if os.path.exists('cpg_matched_probes.txt'):
    with open('cpg_matched_probes.txt', 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:
                probe = parts[0]
                if 'CHI3L1' in line:
                    chi3l1_probes.append(probe)
                elif 'CD44' in line:
                    cd44_probes.append(probe)
print(f"  Probes found in reference: CHI3L1={len(chi3l1_probes)}, CD44={len(cd44_probes)}")

# ── STEP 2: Load Methylation Data ──
print("\nSTEP 2: Loading methylation data for target promoter probes...")
meth = pd.read_csv(os.path.join(data_dir, "beta.450.tsv"), sep='\t', index_col=0)
all_target_probes = chi3l1_probes + cd44_probes
meth_t = meth.loc[meth.index.intersection(all_target_probes)]
del meth; gc.collect()
print(f"  Loaded methylation values for {len(meth_t)} probes across {len(meth_t.columns)} samples.")

# ── STEP 3: Find Matched Methylation Pairs (TP -> R1) ──
print("\nSTEP 3: Finding matched TP -> R1 methylation pairs...")
meth_cols = []
for col in meth_t.columns:
    # Methylation columns use hyphens: TCGA-DH-A669-TP-12M-450-QDG6WF
    p = col.split('-')
    if len(p) >= 4:
        case = '-'.join(p[:3])
        meth_cols.append({'sample': col, 'case': case, 'tp': p[3]})
mdf = pd.DataFrame(meth_cols)

meth_matched = {}
for case in mdf['case'].unique():
    tp_samples = mdf[(mdf['case'] == case) & (mdf['tp'] == 'TP')]['sample'].tolist()
    r1_samples = mdf[(mdf['case'] == case) & (mdf['tp'] == 'R1')]['sample'].tolist()
    if tp_samples and r1_samples:
        meth_matched[case] = {
            'tp': tp_samples[0],
            'r1': r1_samples[0]
        }
print(f"  Matched TP->R1 methylation pairs: {len(meth_matched)}")

# ── STEP 4: Load matched RNA fold changes ──
print("\nSTEP 4: Loading matched RNA-seq data...")
tpm = pd.read_csv(os.path.join(data_dir, "gene_tpm_matrix_all_samples.tsv"), sep='\t', index_col=0)
targets = ['CHI3L1', 'CD44']
tpm_t = tpm.loc[targets]
del tpm; gc.collect()

rna_cols = []
for col in tpm_t.columns:
    p = col.split('.')
    if len(p) >= 4:
        case = '-'.join(p[:3])
        rna_cols.append({'sample': col, 'case': case, 'tp': p[3]})
rdf = pd.DataFrame(rna_cols)

rna_matched = {}
for case in rdf['case'].unique():
    tp_samples = rdf[(rdf['case'] == case) & (rdf['tp'] == 'TP')]['sample'].tolist()
    r1_samples = rdf[(rdf['case'] == case) & (rdf['tp'] == 'R1')]['sample'].tolist()
    if tp_samples and r1_samples:
        rna_matched[case] = {
            'tp': tp_samples[0],
            'r1': r1_samples[0]
        }

# Find common patients with both matched RNA-seq and methylation
common_cases = list(set(meth_matched.keys()).intersection(rna_matched.keys()))
print(f"  Patients with matched TP->R1 RNA AND Methylation: {len(common_cases)}")

# ── STEP 5: Correlate delta Beta with delta RNA ──
print("\nSTEP 5: Testing promoter methylation change vs expression change...")
results = []
for case in common_cases:
    # Calculate RNA log2FC
    tp_rna_s = rna_matched[case]['tp']
    r1_rna_s = rna_matched[case]['r1']
    
    # Calculate methylation delta
    tp_meth_s = meth_matched[case]['tp']
    r1_meth_s = meth_matched[case]['r1']
    
    for g, probes in [('CHI3L1', chi3l1_probes), ('CD44', cd44_probes)]:
        rna_tp = tpm_t.loc[g, tp_rna_s]
        rna_r1 = tpm_t.loc[g, r1_rna_s]
        delta_rna = np.log2((rna_r1 + 1) / (rna_tp + 1))
        
        for p in probes:
            if p in meth_t.index:
                val_tp = meth_t.loc[p, tp_meth_s]
                val_r1 = meth_t.loc[p, r1_meth_s]
                if not np.isnan(val_tp) and not np.isnan(val_r1):
                    delta_beta = val_r1 - val_tp
                    results.append({
                        'case': case,
                        'gene': g,
                        'probe': p,
                        'delta_beta': delta_beta,
                        'delta_rna': delta_rna,
                        'beta_tp': val_tp,
                        'beta_r1': val_r1
                    })

res_df = pd.DataFrame(results)
if len(res_df) > 0:
    res_df.to_csv(os.path.join(out_dir, "stage5_epigenome_expression_raw.txt"), sep='\t', index=False)

    # Run correlations per probe
    summary_stats = []
    for (gene, probe), group in res_df.groupby(['gene', 'probe']):
        if len(group) >= 5:
            r, p = stats.spearmanr(group['delta_beta'], group['delta_rna'])
            summary_stats.append({
                'gene': gene,
                'probe': probe,
                'n_patients': len(group),
                'correlation': round(r, 4),
                'p_value': p,
                'mean_beta_tp': round(group['beta_tp'].mean(), 4),
                'mean_beta_r1': round(group['beta_r1'].mean(), 4),
                'mean_delta_beta': round(group['delta_beta'].mean(), 4),
            })

    summary_df = pd.DataFrame(summary_stats)
    if len(summary_df) > 0:
        summary_df = summary_df.sort_values('p_value')
        summary_df.to_csv(os.path.join(out_dir, "stage5_probe_correlation_summary.txt"), sep='\t', index=False)
        print("\n=== STAGE 5 CORRELATION RESULTS ===")
        for _, r in summary_df.iterrows():
            sig_str = "*" if r['p_value'] < 0.05 else ""
            print(f"  {r['gene']:7s} | {r['probe']} | n={r['n_patients']:2d} | mean_delta_beta={r['mean_delta_beta']:+.3f} | r={r['correlation']:+.3f} | p={r['p_value']:.2e}{sig_str}")
    else:
        print("  No testable probes found.")
else:
    print("  No matched cases with data.")

print("\nStage 5 complete.")
