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

# Load matched variants
df_matches = pd.read_csv(r"D:\research of the GBM\clinvar_pathogenic_matches.txt", sep="\t")

# Load ClinVar summaries to extract GRCh38 coordinates
with open(r"D:\research of the GBM\clinvar_summaries.json") as f:
    clinvar_data = json.load(f)

grch38_lookup = {}
for entry in clinvar_data:
    uid = entry['uid']
    raw = entry['raw_summary']
    var_set = raw.get('variation_set', [])
    if var_set:
        locs = var_set[0].get('variation_loc', [])
        for loc in locs:
            if loc.get('assembly_name') == 'GRCh38':
                grch38_lookup[uid] = {
                    'chrom': str(loc.get('chr')),
                    'pos': int(loc.get('start'))
                }

# Add GRCh38 coordinates to matches
matches_with_coords = []
for idx, row in df_matches.iterrows():
    uid = str(row['clinvar_uid'])
    if uid in grch38_lookup:
        c = grch38_lookup[uid]
        matches_with_coords.append({
            'variant_id': row['variant_id'],
            'gene_symbol': row['gene_symbol'],
            'chrom_grch38': 'chr' + c['chrom'],
            'pos_grch38': c['pos'],
            'ref': row['ref'],
            'alt': row['alt'],
            'clinvar_clinsig': row['clinvar_clinsig']
        })

df_coords = pd.DataFrame(matches_with_coords)
print(f"Loaded {len(df_coords)} variants with GRCh38 coordinates.")

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

# Run a subset of variants first to check (e.g. all 3 IDH1 variants, all 26 ATRX variants, and 20 TP53 variants)
# This keeps it fast and handles rate limits/timeouts
target_variants = pd.concat([
    df_coords[df_coords['gene_symbol'] == 'IDH1'],
    df_coords[df_coords['gene_symbol'] == 'ATRX'],
    df_coords[df_coords['gene_symbol'] == 'TP53'].head(20)
])

print(f"Scoring {len(target_variants)} selected pathogenic variants...")

for idx, row in target_variants.iterrows():
    chrom = row['chrom_grch38']
    pos = row['pos_grch38']
    # If ref/alt are longer than 1 (indels), get first char to avoid API errors
    ref = str(row['ref'])[0]
    alt = str(row['alt'])[0]
    
    interval = genome.Interval(chrom, pos - SEQ_LENGTH // 2, pos + SEQ_LENGTH // 2)
    variant = genome.Variant(chrom, pos, ref, alt)
    
    print(f"  Scoring {row['gene_symbol']} variant at {chrom}:{pos}...")
    try:
        scores_list = dna_model.score_variant(interval=interval, variant=variant, variant_scorers=scorers)
        
        # Process and tidy scores
        max_dnase_quantile = 0.0
        max_histone_quantile = 0.0
        max_tf_quantile = 0.0
        
        for score_adata in scores_list:
            df_tidy = variant_scorers.tidy_scores([score_adata], match_gene_strand=True)
            if df_tidy is not None and not df_tidy.empty:
                # Find max quantile score per modality
                # In some versions it uses 'modality' or 'output_type'
                # Let's inspect unique values of output_type/modality
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
            'variant_id': row['variant_id'],
            'gene_symbol': row['gene_symbol'],
            'chrom': chrom,
            'pos': pos,
            'ref': ref,
            'alt': alt,
            'max_dnase_quantile': max_dnase_quantile,
            'max_histone_quantile': max_histone_quantile,
            'max_tf_quantile': max_tf_quantile
        })
    except Exception as e:
        print(f"    Error scoring variant: {e}")

df_results = pd.DataFrame(results)
df_results.to_csv(r"D:\research of the GBM\alphagenome_pathogenicity_scores.txt", sep="\t", index=False)
print("\nPredictions complete. Summary of results:")
print(df_results.to_string())
