"""CGGA Cross-Validation of Top Probes.
Correlates baseline CpG methylation with matched RNA expression in CGGA IDH-A cohort.
"""
import pandas as pd, numpy as np, os, gc, sys
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8')

cgga_dir = r"D:\research of the GBM"

# Top probes from Stage 5 matched check
target_probes = {
    'cg13134650': 'CHI3L1',
    'cg15427520': 'CD44',
    'cg17014757': 'CHI3L1',
    'cg15490070': 'CHI3L1',
    'cg03625911': 'CHI3L1',
    'cg08530414': 'CD44',
    'cg04361579': 'CHI3L1'
}

# ── STEP 1: Load CGGA Methylation for Target Probes ──
print("STEP 1: Loading CGGA methylation matrix...")
meth = pd.read_csv(os.path.join(cgga_dir, "CGGA_IDH_A_Methylation_EPIC_Array_20250915.txt"), sep='\t', index_col=0)
meth_t = meth.loc[meth.index.intersection(target_probes.keys())]
del meth; gc.collect()
print(f"  Loaded {len(meth_t)} target probes across {len(meth_t.columns)} samples.")

# ── STEP 2: Load CGGA RNA-seq for Target Genes ──
print("\nSTEP 2: Loading CGGA RNA-seq expression...")
rna = pd.read_csv(os.path.join(cgga_dir, "CGGA_IDH_A_RNAseq_RSEM_20250915.txt"), sep='\t', index_col=0)
targets = ['CHI3L1', 'CD44']
rna_t = rna.loc[rna.index.intersection(targets)]
del rna; gc.collect()
print(f"  Loaded expression for {len(rna_t)} target genes across {len(rna_t.columns)} samples.")

# ── STEP 3: Find common patients ──
common_samples = list(set(meth_t.columns).intersection(rna_t.columns))
print(f"\nSTEP 3: Mapped common samples with both data: {len(common_samples)}")
print("Common samples:", common_samples)

# ── STEP 4: Correlate Baseline Beta with Expression ──
print("\nSTEP 4: Testing CpG methylation vs RNA expression correlation...")
results = []
for p, g in target_probes.items():
    if p in meth_t.index and g in rna_t.index:
        meth_vals = meth_t.loc[p, common_samples].values.astype(float)
        rna_vals = rna_t.loc[g, common_samples].values.astype(float)
        
        # Filter NaNs
        valid = ~np.isnan(meth_vals) & ~np.isnan(rna_vals)
        if valid.sum() >= 5:
            r, p_val = stats.spearmanr(meth_vals[valid], rna_vals[valid])
            results.append({
                'gene': g,
                'probe': p,
                'n_samples': int(valid.sum()),
                'correlation': round(r, 4),
                'p_value': p_val
            })

df_res = pd.DataFrame(results)
if len(df_res) > 0:
    df_res = df_res.sort_values('p_value')
    df_res.to_csv(os.path.join(cgga_dir, "cgga_cpg_expression_cross_val.txt"), sep='\t', index=False)
    print("\n=== CGGA CROSS-VALIDATION RESULTS ===")
    for _, r in df_res.iterrows():
        sig_str = "*" if r['p_value'] < 0.05 else ""
        print(f"  {r['gene']:7s} | {r['probe']} | n={r['n_samples']:2d} | r={r['correlation']:+.3f} | p={r['p_value']:.2e}{sig_str}")
else:
    print("  No correlations could be calculated.")
