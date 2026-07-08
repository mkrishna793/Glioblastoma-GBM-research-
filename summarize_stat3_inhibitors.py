import json
import pandas as pd

file_path = r"C:\Users\bhanu\.gemini\antigravity\scratch\chembl_stat3_activities.json"
data = json.load(open(file_path))

acts = data.get('activities', [])
records = []
for a in acts:
    mol_id = a.get('molecule_chembl_id')
    std_val = a.get('standard_value')
    std_type = a.get('standard_type')
    std_units = a.get('standard_units')
    desc = a.get('assay_description', '')
    
    if std_val is not None:
        try:
            val = float(std_val)
            records.append({
                'Molecule_ID': mol_id,
                'IC50_nM': val,
                'Assay': desc[:120] + '...' if len(desc) > 120 else desc
            })
        except ValueError:
            pass

df = pd.DataFrame(records)
if not df.empty:
    df_sorted = df.sort_values(by='IC50_nM', ascending=True)
    print("Top 10 most potent STAT3 inhibitors in ChEMBL:")
    print(df_sorted.head(10).to_string(index=False))
else:
    print("No valid bioactivity records found.")
