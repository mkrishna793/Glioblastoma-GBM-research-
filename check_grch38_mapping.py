import json
import pandas as pd

print("Loading matched variants...")
df_matches = pd.read_csv(r"D:\research of the GBM\clinvar_pathogenic_matches.txt", sep="\t")

print("Loading ClinVar summaries...")
with open(r"D:\research of the GBM\clinvar_summaries.json") as f:
    clinvar_data = json.load(f)

# Build a lookup for GRCh38 coordinates by ClinVar Variation ID
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
                    'pos': loc.get('start'),
                    # Reference and alternate alleles are not directly in variation_loc in some summaries,
                    # but we can look for them or use the local ones (or HGVS)
                }

# Check how many matched variants have GRCh38 coordinates
count = 0
for idx, row in df_matches.iterrows():
    uid = str(row['clinvar_uid'])
    if uid in grch38_lookup:
        count += 1
print(f"Matched variants with GRCh38 coordinates: {count} / {len(df_matches)}")
