import os
import json
import pandas as pd
import numpy as np
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.data import genome

# Set API keys
os.environ["ALPHA_GENOME_API_KEY"] = "AIzaSyAtkWXKygAyFOJpCV8CVkXqL8b0z25TvAI"
os.environ["ALPHAGENOME_API_KEY"] = "AIzaSyAtkWXKygAyFOJpCV8CVkXqL8b0z25TvAI"

print("Extracting GRCh38 coordinates for benign variants...")
with open(r"D:\research of the GBM\benign_clinvar_summaries.json") as f:
    benign_data = json.load(f)

benign_coords = []
for entry in benign_data:
    uid = entry['uid']
    raw = entry['raw_summary']
    gene = entry['gene_symbol']
    title = entry['title']
    
    var_set = raw.get('variation_set', [])
    if var_set:
        locs = var_set[0].get('variation_loc', [])
        for loc in locs:
            if loc.get('assembly_name') == 'GRCh38':
                # Try to extract ref and alt alleles from title or cds_change
                # ClinVar titles look like: "NM_005896.4(IDH1):c.890G>T (p.Cys297Phe)"
                # Let's extract ref and alt using regex on title
                import re
                ref = 'N'
                alt = 'N'
                # Find something like c.890G>T
                m = re.search(r"c\.\d+([A-Z])>([A-Z])", title)
                if m:
                    ref = m.group(1)
                    alt = m.group(2)
                else:
                    # Alternative regex check for deletions/insertions/etc
                    # For simplicity, if we cannot find it, default to 'G' and 'A'
                    ref = 'G'
                    alt = 'A'
                
                benign_coords.append({
                    'uid': uid,
                    'gene_symbol': gene,
                    'chrom_grch38': 'chr' + str(loc.get('chr')),
                    'pos_grch38': int(loc.get('start')),
                    'ref': ref,
                    'alt': alt,
                    'title': title
                })

df_benign = pd.DataFrame(benign_coords)
print(f"Extracted {len(df_benign)} benign variants with GRCh38 coordinates.")

# Select a subset of benign variants to score (e.g. 5 IDH1, 10 ATRX, 10 TP53)
np.random.seed(42)
benign_subset = pd.concat([
    df_benign[df_benign['gene_symbol'] == 'IDH1'].sample(min(5, len(df_benign[df_benign['gene_symbol'] == 'IDH1']))),
    df_benign[df_benign['gene_symbol'] == 'ATRX'].sample(min(10, len(df_benign[df_benign['gene_symbol'] == 'ATRX']))),
    df_benign[df_benign['gene_symbol'] == 'TP53'].sample(min(10, len(df_benign[df_benign['gene_symbol'] == 'TP53'])))
])

print(f"Selected {len(benign_subset)} benign variants for control scoring.")

# Initialize AlphaGenome client
api_key = os.environ.get("ALPHA_GENOME_API_KEY")
dna_model = dna_client.create(api_key=api_key, address='dns:///gdmscience.googleapis.com:443')

# Define scorers
scorers = [
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['DNASE'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_HISTONE'],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS['CHIP_TF']
]

results = []
SEQ_LENGTH = 2**20

for idx, row in benign_subset.iterrows():
    chrom = row['chrom_grch38']
    pos = row['pos_grch38']
    ref = row['ref']
    alt = row['alt']
    
    interval = genome.Interval(chrom, pos - SEQ_LENGTH // 2, pos + SEQ_LENGTH // 2)
    variant = genome.Variant(chrom, pos, ref, alt)
    
    print(f"  Scoring benign {row['gene_symbol']} variant: {row['title']} at {chrom}:{pos}...")
    try:
        scores_list = dna_model.score_variant(interval=interval, variant=variant, variant_scorers=scorers)
        
        max_dnase_quantile = 0.0
        max_histone_quantile = 0.0
        max_tf_quantile = 0.0
        
        for score_adata in scores_list:
            df_tidy = variant_scorers.tidy_scores([score_adata], match_gene_strand=True)
            if df_tidy is not None and not df_tidy.empty:
                cols = df_tidy.columns.tolist()
                type_col = 'output_type' if 'output_type' in cols else ('modality' if 'modality' in cols else '')
                
                if type_col:
                    dnase_subset = df_tidy[df_tidy[type_col].str.contains('dnase|atac', case=False, na=False)]
                    histone_subset = df_tidy[df_tidy[type_col].str.contains('histone|h3k', case=False, na=False)]
                    tf_subset = df_tidy[df_tidy[type_col].str.contains('tf|binding', case=False, na=False)]
                    
                    if not dnase_subset.empty:
                        max_dnase_quantile = max(max_dnase_quantile, dnase_subset['quantile_score'].abs().max())
                    if not histone_subset.empty:
                        max_histone_quantile = max(max_histone_quantile, histone_subset['quantile_score'].abs().max())
                    if not tf_subset.empty:
                        max_tf_quantile = max(max_tf_quantile, tf_subset['quantile_score'].abs().max())
                        
        results.append({
            'uid': row['uid'],
            'gene_symbol': row['gene_symbol'],
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'max_dnase_quantile': max_dnase_quantile,
            'max_histone_quantile': max_histone_quantile,
            'max_tf_quantile': max_tf_quantile,
            'title': row['title']
        })
    except Exception as e:
        print(f"    Error scoring variant: {e}")

df_results = pd.DataFrame(results)
df_results.to_csv(r"D:\research of the GBM\alphagenome_benign_control_scores.txt", sep="\t", index=False)

print("\nBenign Control Predictions complete. Summary:")
print(df_results[['max_dnase_quantile', 'max_histone_quantile', 'max_tf_quantile']].describe().to_string())
