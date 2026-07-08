import urllib.request
import json
import time

# Load benign IDs
with open(r"D:\research of the GBM\benign_clinvar_ids.json") as f:
    benign_ids = json.load(f)

all_ids = []
for g in benign_ids:
    all_ids.extend(benign_ids[g])
print(f"Total benign IDs to fetch: {len(all_ids)}")

# Fetch summaries
benign_summaries = []
batch_size = 100
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
                title = summary.get('title', '')
                genes_list = summary.get('genes', [])
                gene_symbol = genes_list[0].get('symbol', '') if genes_list else ''
                
                benign_summaries.append({
                    'uid': uid,
                    'title': title,
                    'gene_symbol': gene_symbol,
                    'raw_summary': summary
                })
        time.sleep(1)
    except Exception as e:
        print(f"  Error fetching batch starting at {i}: {e}")

with open(r"D:\research of the GBM\benign_clinvar_summaries.json", "w") as f:
    json.dump(benign_summaries, f, indent=2)
print("Saved benign summaries to D:\\research of the GBM\\benign_clinvar_summaries.json")
