import json
import pandas as pd
import re

print("Loading local variants...")
df_local = pd.read_csv(r"D:\research of the GBM\key_genes_coordinates.txt", sep="\t")

print("Loading ClinVar summaries...")
with open(r"D:\research of the GBM\clinvar_summaries.json") as f:
    clinvar_data = json.load(f)

# Build a lookup dictionary of ClinVar variants by GRCh37 coordinate
# Key: (gene_symbol, chrom, position_1indexed) -> ClinVar record details
clinvar_lookup = {}
for entry in clinvar_data:
    raw = entry['raw_summary']
    gene = entry['gene_symbol']
    title = entry['title']
    clinsig = raw.get('clinical_significance', {}).get('description', '')
    
    var_set = raw.get('variation_set', [])
    if var_set:
        locs = var_set[0].get('variation_loc', [])
        for loc in locs:
            if loc.get('assembly_name') == 'GRCh37':
                chr_num = str(loc.get('chr'))
                start = loc.get('start')
                if start is not None:
                    try:
                        pos = int(start)
                        key = (gene, chr_num, pos)
                        clinvar_lookup[key] = {
                            'title': title,
                            'clinsig': clinsig,
                            'uid': entry['uid']
                        }
                    except ValueError:
                        pass

print(f"Built ClinVar lookup index with {len(clinvar_lookup)} GRCh37 variants.")

# Match each local variant
matches = []
for idx, row in df_local.iterrows():
    # Parse pos range e.g. [209103840,209103841)
    pos_str = row['pos']
    m = re.match(r"\[(\d+),(\d+)\)", pos_str)
    if m:
        # Range is [start, end). The 1-indexed position is start + 1
        pos_1indexed = int(m.group(1)) + 1
        gene = row['gene_symbol']
        chrom = str(row['chrom'])
        if chrom == '23':
            chrom = 'X'
        elif chrom == '24':
            chrom = 'Y'
        
        # Check lookup
        key = (gene, chrom, pos_1indexed)
        if key in clinvar_lookup:
            matches.append({
                'variant_id': row['variant_id'],
                'gene_symbol': gene,
                'chrom': chrom,
                'pos_local': pos_str,
                'pos_1indexed': pos_1indexed,
                'ref': row['ref'],
                'alt': row['alt'],
                'variant_classification': row['variant_classification'],
                'clinvar_title': clinvar_lookup[key]['title'],
                'clinvar_clinsig': clinvar_lookup[key]['clinsig'],
                'clinvar_uid': clinvar_lookup[key]['uid']
            })

df_matches = pd.DataFrame(matches)
print(f"Matched {len(df_matches)} local variants to ClinVar pathogenic variants!")
if not df_matches.empty:
    print(df_matches[['gene_symbol', 'pos_1indexed', 'ref', 'alt', 'clinvar_clinsig']].value_counts())
    df_matches.to_csv(r"D:\research of the GBM\clinvar_pathogenic_matches.txt", sep="\t", index=False)
