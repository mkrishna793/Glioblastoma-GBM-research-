import subprocess
import json
import os

chembl_script = r"C:\Users\bhanu\.gemini\config\plugins\science\skills\chembl_database\scripts\chembl_api.py"
output_dir = r"D:\research of the GBM"

# Helper to run chembl_api.py
def run_chembl_cmd(subcommand, args_dict):
    cmd = ["uv", "run", chembl_script, subcommand]
    for k, v in args_dict.items():
        cmd.extend([k, str(v)])
    
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        return None
    return result.stdout

# 1. Search for targets: CHI3L1, STAT3, NFKB1, IKBKB, TEAD1, YAP1
targets_to_search = ["CHI3L1", "STAT3", "NFKB1", "IKBKB", "TEAD1", "YAP1"]
target_results = {}

for t in targets_to_search:
    out_file = os.path.join(output_dir, f"chembl_target_{t}.json")
    args = {
        "--search": t,
        "--limit": 5,
        "--output": out_file
    }
    run_chembl_cmd("target", args)
    
    if os.path.exists(out_file):
        with open(out_file, 'r') as f:
            try:
                data = json.load(f)
                target_results[t] = data.get("targets", [])
                print(f"Found {len(target_results[t])} targets for {t}")
            except Exception as e:
                print(f"Error reading {out_file}: {e}")

# Save summarized targets
with open(os.path.join(output_dir, "chembl_target_summaries.json"), "w") as f:
    json.dump(target_results, f, indent=2)

# 2. Get activities / inhibitors for key targets
# Let's find activities for STAT3 (CHEMBL5147) and IKBKB (CHEMBL2007) and CHI3L1 if exists
# We will check target IDs from the search results
for t, targets in target_results.items():
    if not targets:
        continue
    # Get the first target's chembl_id
    target_id = targets[0]["target_chembl_id"]
    target_name = targets[0].get("pref_name", t)
    print(f"Fetching activities for target {target_name} ({target_id})...")
    
    out_act_file = os.path.join(output_dir, f"chembl_activity_{t}.json")
    act_args = {
        "--filter": f"target_chembl_id={target_id} standard_type=IC50",
        "--limit": 10,
        "--normalize": "",
        "--output": out_act_file
    }
    run_chembl_cmd("activity", act_args)

print("ChEMBL queries complete.")
