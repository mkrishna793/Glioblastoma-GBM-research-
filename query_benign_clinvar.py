import urllib.request
import json
import time

print("Querying ClinVar for Benign/Likely Benign variants in IDH1, ATRX, TP53...")
genes = ['IDH1', 'ATRX', 'TP53']
benign_ids = {}

for gene in genes:
    term = f"{gene}[gene] AND (benign[clinsig] OR likely benign[clinsig])"
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={urllib.parse.quote(term)}&retmode=json&retmax=100" # We only need a representative set of 100 max per gene
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            res = json.loads(response.read().decode('utf-8'))
            ids = res.get('esearchresult', {}).get('idlist', [])
            benign_ids[gene] = ids
            print(f"  {gene}: found {len(ids)} benign/likely benign variant IDs in ClinVar")
            time.sleep(1) # Polite pause
    except Exception as e:
        print(f"  Error querying {gene}: {e}")

# Save the IDs
with open(r"D:\research of the GBM\benign_clinvar_ids.json", "w") as f:
    json.dump(benign_ids, f, indent=2)
print("Saved benign IDs to D:\\research of the GBM\\benign_clinvar_ids.json")
