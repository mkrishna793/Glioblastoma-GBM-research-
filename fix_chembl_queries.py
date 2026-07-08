import subprocess
import json
import os

chembl_script = r"C:\Users\bhanu\.gemini\config\plugins\science\skills\chembl_database\scripts\chembl_api.py"
output_dir = r"D:\research of the GBM"

# Targets to search
targets = ["CHI3L1", "STAT3", "NFKB1", "IKBKB", "TEAD1", "YAP1"]

for t in targets:
    in_file = os.path.join(output_dir, f"chembl_target_{t}.json")
    if not os.path.exists(in_file):
        continue
    
    with open(in_file, 'r') as f:
        try:
            data = json.load(f)
            target_list = data.get("targets", [])
            if not target_list:
                print(f"No targets found for {t}")
                continue
            
            # Find the best target (human target)
            best_target = None
            for tgt in target_list:
                if tgt.get("organism") == "Homo sapiens" and tgt.get("target_type") == "SINGLE PROTEIN":
                    best_target = tgt
                    break
            
            if not best_target:
                best_target = target_list[0]
            
            target_id = best_target["target_chembl_id"]
            pref_name = best_target.get("pref_name", t)
            print(f"Target found for {t}: {pref_name} ({target_id})")
            
            # Query activities using correct arguments list
            out_act_file = os.path.join(output_dir, f"chembl_activity_fixed_{t}.json")
            cmd = [
                "uv", "run", chembl_script, "activity",
                "--filter", f"target_chembl_id={target_id} standard_type=IC50",
                "--limit", "10",
                "--normalize",
                "--output", out_act_file
            ]
            
            print(f"Running command: {' '.join(cmd)}")
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                print(f"Error querying activity for {t}: {res.stderr}")
            else:
                print(f"Successfully saved activity file for {t}")
                
        except Exception as e:
            print(f"Error parsing {in_file}: {e}")
