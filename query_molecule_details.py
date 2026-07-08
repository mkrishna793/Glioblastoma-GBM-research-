import subprocess
import json
import os

chembl_script = r"C:\Users\bhanu\.gemini\config\plugins\science\skills\chembl_database\scripts\chembl_api.py"
output_dir = r"D:\research of the GBM"

mols = ["CHEMBL5568901", "CHEMBL5401444", "CHEMBL502473", "CHEMBL4439905"]
out_file = os.path.join(output_dir, "chembl_molecule_details.json")

# Semicolon-separated string for batch fetch
ids_str = ";".join(mols)

cmd = [
    "uv", "run", chembl_script, "molecule",
    "--ids", ids_str,
    "--limit", "10",
    "--output", out_file
]

print(f"Running: {' '.join(cmd)}")
res = subprocess.run(cmd, capture_output=True, text=True)
if res.returncode != 0:
    print(f"Error: {res.stderr}")
else:
    if os.path.exists(out_file):
        with open(out_file, 'r') as f:
            try:
                data = json.load(f)
                molecules = data.get("molecules", [])
                print(f"Successfully fetched details for {len(molecules)} molecules:")
                for m in molecules:
                    cid = m.get("molecule_chembl_id")
                    pref_name = m.get("pref_name", "None")
                    phase = m.get("max_phase", "0")
                    syns = [s.get("molecule_synonym", "") for s in m.get("molecule_synonyms", []) if s.get("syn_type") == "TRADE_NAME" or s.get("syn_type") == "USP"]
                    print(f"  ID: {cid} | Name: {pref_name} | Max Clinical Phase: {phase} | Synonyms: {', '.join(syns[:2])}")
            except Exception as e:
                print(f"Error reading molecule file: {e}")
