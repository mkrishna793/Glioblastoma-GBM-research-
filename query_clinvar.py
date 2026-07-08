import urllib.request
import json
import time

print("Querying ClinVar via NCBI E-utilities (urllib)...")
genes = ['IDH1', 'ATRX', 'TP53']
pathogenic_ids = {}

for gene in genes:
    term = f"{gene}[gene] AND (pathogenic[clinsig] OR likely pathogenic[clinsig])"
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={urllib.parse.quote(term)}&retmode=json&retmax=1000"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            res = json.loads(response.read().decode('utf-8'))
            ids = res.get('esearchresult', {}).get('idlist', [])
            pathogenic_ids[gene] = ids
            print(f"  {gene}: found {len(ids)} pathogenic/likely pathogenic variant IDs in ClinVar")
            time.sleep(1) # Polite pause
    except Exception as e:
        print(f"  Error querying {gene}: {e}")

# Save the IDs
with open(r"D:\research of the GBM\pathogenic_clinvar_ids.json", "w") as f:
    json.dump(pathogenic_ids, f, indent=2)
print("Saved IDs to D:\\research of the GBM\\pathogenic_clinvar_ids.json")
