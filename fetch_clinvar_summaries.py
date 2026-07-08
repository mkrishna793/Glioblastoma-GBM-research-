import urllib.request
import json
import time
import pandas as pd
import re

print("Loading local variants...")
df_local = pd.read_csv(r"D:\research of the GBM\key_genes_coordinates.txt", sep="\t")
print(f"Loaded {len(df_local)} local variants.")

# Load pathogenic ClinVar IDs
with open(r"D:\research of the GBM\pathogenic_clinvar_ids.json") as f:
    pathogenic_ids = json.load(f)

# Combine all IDs to fetch
all_ids = []
for g in pathogenic_ids:
    all_ids.extend(pathogenic_ids[g])
print(f"Total ClinVar IDs to fetch: {len(all_ids)}")

# Fetch summaries in batches of 300
clinvar_coords = []
batch_size = 300
for i in range(0, len(all_ids), batch_size):
    batch = all_ids[i:i+batch_size]
    ids_str = ",".join(batch)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={ids_str}&retmode=json"
    
    print(f"  Fetching batch {i//batch_size + 1}/{len(all_ids)//batch_size + 1}...")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            res = json.loads(response.read().decode('utf-8'))
            result = res.get('result', {})
            for uid in batch:
                summary = result.get(uid, {})
                # Extract coordinates
                # ClinVar JSON summary contains 'variation_loc' or similar, let's extract it
                # We can print one summary to inspect structure if needed
                title = summary.get('title', '')
                genes_list = summary.get('genes', [])
                gene_symbol = genes_list[0].get('symbol', '') if genes_list else ''
                
                # Check for genomic coordinates in the summary
                # Let's collect the raw summary to parse it in memory
                clinvar_coords.append({
                    'uid': uid,
                    'title': title,
                    'gene_symbol': gene_symbol,
                    'raw_summary': summary
                })
        time.sleep(1)
    except Exception as e:
        print(f"  Error fetching batch starting at {i}: {e}")

# Save raw ClinVar summaries
with open(r"D:\research of the GBM\clinvar_summaries.json", "w") as f:
    json.dump(clinvar_coords, f, indent=2)
print("Saved ClinVar summaries.")
