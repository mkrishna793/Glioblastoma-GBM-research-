"""Stage 1: Unbiased Mutation Scan — genes driving CHI3L1/CD44 at recurrence."""
import pandas as pd, numpy as np, os, gc, sys
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8')

data_dir = r"D:\new data of the GBM research"
out_dir  = r"D:\research of the GBM"

def fdr_bh(pvals):
    n = len(pvals)
    idx = np.argsort(pvals)
    adj = np.array(pvals)[idx] * n / np.arange(1, n+1)
    for i in range(n-2, -1, -1):
        adj[i] = min(adj[i], adj[i+1])
    result = np.empty(n)
    result[idx] = np.minimum(adj, 1.0)
    return result

# ── STEP 1: Build coding variant lookup ──
print("STEP 1: Loading coding variant annotations...")
coding_cls = {
    'MISSENSE', 'NONSENSE', 'FRAME_SHIFT_DEL', 'FRAME_SHIFT_INS',
    'SPLICE_SITE', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'NONSTOP',
    'DE_NOVO_START_IN_FRAME', 'DE_NOVO_START_OUT_FRAME',
    'START_CODON_SNP', 'START_CODON_DEL', 'START_CODON_INS'
}

vid_gene = {}
for chunk in pd.read_csv(os.path.join(data_dir, "variants.anno.csv.gz"),
        usecols=['variant_id','gene_symbol','variant_classification'], chunksize=500000):
    c = chunk[chunk['variant_classification'].isin(coding_cls)]
    vid_gene.update(dict(zip(c['variant_id'], c['gene_symbol'])))
print(f"  Coding variant IDs loaded: {len(vid_gene)}")
gc.collect()

# ── STEP 2: Build patient × gene mutation matrix ──
print("\nSTEP 2: Building patient-gene matrix...")
patient_genes = {}
patient_tp = {}
artifact_re = r'\.\d+$|_pyclone'
n_hits = 0

for chunk in pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"),
        usecols=['aliquot_barcode','variant_id','case_barcode'], chunksize=300000):
    chunk = chunk[~chunk['aliquot_barcode'].astype(str).str.contains(artifact_re, regex=True, na=False)]
    chunk['gene'] = chunk['variant_id'].map(vid_gene)
    chunk = chunk.dropna(subset=['gene'])
    n_hits += len(chunk)
    chunk['tp'] = chunk['aliquot_barcode'].str.split('-').str[3]
    for case, gene, tp in zip(chunk['case_barcode'], chunk['gene'], chunk['tp']):
        patient_genes.setdefault(case, set()).add(gene)
        patient_tp.setdefault(case, set()).add(tp)

print(f"  Coding genotype hits: {n_hits}")
print(f"  Patients with coding muts: {len(patient_genes)}")
del vid_gene; gc.collect()

# ── STEP 3: RNA fold-change TP->R1 for CHI3L1/CD44 ──
print("\nSTEP 3: Computing RNA fold-change TP->R1...")
tpm = pd.read_csv(os.path.join(data_dir, "gene_tpm_matrix_all_samples.tsv"), sep='\t', index_col=0)
targets = [g for g in ['CHI3L1','CD44'] if g in tpm.index]
tpm_t = tpm.loc[targets]
del tpm; gc.collect()

cols_info = []
for col in tpm_t.columns:
    # TPM columns use dots (e.g. GLSS.19.0266.R1.01R.RNA.FFKUM5)
    p = col.split('.')
    if len(p) >= 4:
        case_hyphenated = '-'.join(p[:3])
        cols_info.append({'sample': col, 'case': case_hyphenated, 'tp': p[3]})
cdf = pd.DataFrame(cols_info)

fold_changes = {}
if len(cdf) > 0:
    for case in cdf['case'].unique():
        tp_s = cdf[(cdf['case']==case)&(cdf['tp']=='TP')]['sample'].tolist()
        r1_s = cdf[(cdf['case']==case)&(cdf['tp']=='R1')]['sample'].tolist()
        if tp_s and r1_s:
            fc = {}
            for g in targets:
                tv = tpm_t.loc[g, tp_s].mean()
                rv = tpm_t.loc[g, r1_s].mean()
                fc[g] = np.log2((rv+1)/(tv+1))
            fold_changes[case] = fc
print(f"  Patients with matched TP->R1 RNA: {len(fold_changes)}")
del tpm_t; gc.collect()

# ── STEP 4: Mann-Whitney U test per gene, FDR correction ──
print("\nSTEP 4: Running gene scan with FDR...")
rna_cases = set(fold_changes.keys())

from collections import Counter
gf = Counter()
for c in rna_cases:
    if c in patient_genes:
        for g in patient_genes[c]:
            gf[g] += 1
testable = [g for g,n in gf.items() if n >= 5]
print(f"  Genes mutated in >=5 patients (with RNA): {len(testable)}")

results = []
for gene in testable:
    mut = {c for c in rna_cases if c in patient_genes and gene in patient_genes[c]}
    wt  = rna_cases - mut
    if len(mut) < 5 or len(wt) < 5:
        continue
    for tgt in targets:
        m_fc = [fold_changes[c][tgt] for c in mut]
        w_fc = [fold_changes[c][tgt] for c in wt]
        U, p = stats.mannwhitneyu(m_fc, w_fc, alternative='two-sided')
        results.append({
            'gene': gene, 'target': tgt,
            'n_mut': len(m_fc), 'n_wt': len(w_fc),
            'med_fc_mut': round(np.median(m_fc),4),
            'med_fc_wt':  round(np.median(w_fc),4),
            'effect':     round(np.median(m_fc)-np.median(w_fc),4),
            'p_value':    p
        })

df = pd.DataFrame(results)
if len(df) > 0:
    df['q_value'] = fdr_bh(df['p_value'].values)
    df['sig'] = df['q_value'] < 0.05
    df = df.sort_values('q_value')
    df.to_csv(os.path.join(out_dir, "stage1_unbiased_scan_all.txt"), sep='\t', index=False)

    sig = df[df['sig']]
    print(f"\n  Total tests: {len(df)}")
    print(f"  Significant after FDR (q<0.05): {len(sig)}")
    if len(sig) > 0:
        print("\n  === FDR-SIGNIFICANT GENES ===")
        for _, r in sig.head(30).iterrows():
            print(f"  {r['gene']:15s} -> {r['target']:7s}  n_mut={r['n_mut']:3d}  effect={r['effect']:+.3f}  p={r['p_value']:.2e}  q={r['q_value']:.4f}")
    print("\n  TOP 20 (regardless of significance):")
    for _, r in df.head(20).iterrows():
        print(f"  {r['gene']:15s} -> {r['target']:7s}  n_mut={r['n_mut']:3d}  effect={r['effect']:+.3f}  p={r['p_value']:.2e}  q={r['q_value']:.4f}")
else:
    print("  No results.")

print("\nStage 1 complete.")
