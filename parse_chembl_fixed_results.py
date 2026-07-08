import json
import os

output_dir = r"D:\research of the GBM"
targets = ["CHI3L1", "STAT3", "NFKB1", "IKBKB", "TEAD1", "YAP1"]

print("=" * 60)
print("EXTRACTING TOP INHIBITORS FROM CHEMBL (FIXED FILTERS)")
print("=" * 60)

for t in targets:
    act_file = os.path.join(output_dir, f"chembl_activity_fixed_filter_{t}.json")
    if not os.path.exists(act_file):
        print(f"No activity file for {t}")
        continue
        
    with open(act_file, 'r') as f:
        try:
            data = json.load(f)
            activities = data.get("activities", [])
            print(f"\n--- Top Compounds for {t} (IC50) ---")
            if not activities:
                print("  No activity records found.")
                continue
                
            # Sort activities by normalized IC50 value (ascending)
            valid_acts = []
            for act in activities:
                val = act.get("normalized_value_nM")
                if val is not None:
                    valid_acts.append(act)
                    
            valid_acts.sort(key=lambda x: float(x["normalized_value_nM"]))
            
            for i, act in enumerate(valid_acts[:5]):
                mol_id = act.get("molecule_chembl_id")
                ic50 = act.get("normalized_value_nM")
                relation = act.get("relation", "=")
                units = "nM"
                print(f"  {i+1}. Molecule: {mol_id} | IC50: {relation} {ic50:.2f} {units}")
                
        except Exception as e:
            print(f"Error parsing {act_file}: {e}")
