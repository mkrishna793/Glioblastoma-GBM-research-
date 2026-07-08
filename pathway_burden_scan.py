"""Stage 1 Pathway-Level Burden Test — testing pathway mutation burden vs CHI3L1/CD44 changes."""
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

# ── STEP 1: Load mappings ──
print("STEP 1: Loading Gene Symbol to Ensembl mapping...")
sym_to_ens = pd.read_csv('symbol_to_ensembl.csv')
sym_to_ens_dict = dict(zip(sym_to_ens['symbol'], sym_to_ens['ensembl_id']))

print("Loading Ensembl to Reactome pathway mapping...")
ens_to_path = {}
with open('Ensembl2Reactome.txt', 'r', encoding='utf-8') as f:
    for line in f:
        if 'Homo sapiens' in line:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                ens_id = parts[0]
                path_id = parts[1]
                path_name = parts[3]
                ens_to_path.setdefault(ens_id, set()).add((path_id, path_name))

print("Loading coding variant annotations...")
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

# ── STEP 2: Build patient × pathway mutation burden matrix ──
print("\nSTEP 2: Mapping patient mutations to Reactome pathways...")
patient_pathways = {} # case_barcode -> pathway_id -> count of mutations
artifact_re = r'\.\d+$|_pyclone'

for chunk in pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"),
        usecols=['aliquot_barcode','variant_id','case_barcode'], chunksize=300000):
    chunk = chunk[~chunk['aliquot_barcode'].astype(str).str.contains(artifact_re, regex=True, na=False)]
    chunk['gene'] = chunk['variant_id'].map(vid_gene)
    chunk = chunk.dropna(subset=['gene'])
    
    # Map to Ensembl, then to pathways
    chunk['ensembl'] = chunk['gene'].map(sym_to_ens_dict)
    chunk = chunk.dropna(subset=['ensembl'])
    
    for case, ens in zip(chunk['case_barcode'], chunk['ensembl']):
        if ens in ens_to_path:
            patient_pathways.setdefault(case, {})
            for pid, pname in ens_to_path[ens]:
                patient_pathways[case][pid] = patient_pathways[case].get(pid, 0) + 1

print(f"  Mapped pathways for {len(patient_pathways)} patients.")
del vid_gene; gc.collect()

# ── STEP 3: Load RNA TP->R1 fold changes ──
print("\nSTEP 3: Loading RNA fold-changes...")
tpm = pd.read_csv(os.path.join(data_dir, "gene_tpm_matrix_all_samples.tsv"), sep='\t', index_col=0)
targets = [g for g in ['CHI3L1','CD44'] if g in tpm.index]
tpm_t = tpm.loc[targets]
del tpm; gc.collect()

cols_info = []
for col in tpm_t.columns:
    p = col.split('.')
    if len(p) >= 4:
        case_hyphenated = '-'.join(p[:3])
        cols_info.append({'sample': col, 'case': case_hyphenated, 'tp': p[3]})
cdf = pd.DataFrame(cols_info)

fold_changes = {}
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

# ── STEP 4: Run Pathway burden scan ──
print("\nSTEP 4: Testing pathway burden association...")
rna_cases = list(fold_changes.keys())

# Get all pathways mutated in at least 5 RNA-matched patients
all_pathways = {}
pathway_names = {}
with open('Ensembl2Reactome.txt', 'r', encoding='utf-8') as f:
    for line in f:
        if 'Homo sapiens' in line:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                pathway_names[parts[1]] = parts[3]

pathway_counts = {}
for c in rna_cases:
    if c in patient_pathways:
        for pid in patient_pathways[c]:
            pathway_counts[pid] = pathway_counts.get(pid, 0) + 1
testable_pathways = [pid for pid, n in pathway_counts.items() if n >= 5]
print(f"  Pathways mutated in >=5 patients: {len(testable_pathways)}")

results = []
for pid in testable_pathways:
    mut_cases = {c for c in rna_cases if c in patient_pathways and pid in patient_pathways[c]}
    wt_cases = set(rna_cases) - mut_cases
    if len(mut_cases) < 5 or len(wt_cases) < 5:
        continue
        
    pname = pathway_names.get(pid, "Unknown")
    for tgt in targets:
        m_fc = [fold_changes[c][tgt] for c in mut_cases]
        w_fc = [fold_changes[c][tgt] for c in wt_cases]
        U, p = stats.mannwhitneyu(m_fc, w_fc, alternative='two-sided')
        results.append({
            'pathway_id': pid, 'pathway_name': pname, 'target': tgt,
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
    df.to_csv(os.path.join(out_dir, "pathway_burden_results_all.txt"), sep='\t', index=False)
    
    sig = df[df['sig']]
    print(f"\n  Total pathways tested: {len(df)}")
    print(f"  Significant after FDR (q<0.05): {len(sig)}")
    
    print("\n  TOP 20 PATHWAYS BY SIGNIFICANCE:")
    for _, r in df.head(20).iterrows():
        print(f"  {r['pathway_name'][:40]:40s} ({r['pathway_id']}) -> {r['target']:7s}  n_mut={r['n_mut']:3d}  effect={r['effect']:+.3f}  p={r['p_value']:.2e}  q={r['q_value']:.4f}")
else:
    print("  No pathway results.")

print("\nPathway burden test complete.")
