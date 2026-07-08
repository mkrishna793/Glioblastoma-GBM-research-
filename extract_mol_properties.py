import json
import os

output_dir = r"D:\research of the GBM"
mol_file = os.path.join(output_dir, "chembl_molecule_details.json")

if os.path.exists(mol_file):
    with open(mol_file, 'r') as f:
        data = json.load(f)
        molecules = data.get("molecules", [])
        print("Detailed Compound Profiles:")
        for m in molecules:
            cid = m.get("molecule_chembl_id")
            props = m.get("molecule_properties", {})
            mw = props.get("mw_freebase", "Unknown")
            alogp = props.get("alogp", "Unknown")
            smiles = m.get("molecule_structures", {}).get("canonical_smiles", "Unknown")
            
            # Identify target category
            target_map = {
                "CHEMBL5568901": "CHI3L1 (YKL-40) Inhibitor",
                "CHEMBL5401444": "TEAD1 (YAP/TAZ co-factor) Inhibitor",
                "CHEMBL502473": "STAT3 Dimerization Inhibitor",
                "CHEMBL4439905": "YAP1-TEAD Interaction Inhibitor"
            }
            category = target_map.get(cid, "Target Unknown")
            
            print(f"\nTarget: {category} ({cid})")
            print(f"  Molecular Weight: {mw} g/mol")
            print(f"  LogP (lipophilicity): {alogp}")
            print(f"  Canonical SMILES: {smiles}")
