import subprocess
import json
import os

chembl_script = r"C:\Users\bhanu\.gemini\config\plugins\science\skills\chembl_database\scripts\chembl_api.py"
output_dir = r"D:\research of the GBM"

target_ids = {
    "CHI3L1": "CHEMBL5724768",
    "STAT3": "CHEMBL4026",
    "NFKB1": "CHEMBL3251",
    "IKBKB": "CHEMBL1991",
    "TEAD1": "CHEMBL3334416",
    "YAP1": "CHEMBL3334415"
}

for name, tid in target_ids.items():
    out_file = os.path.join(output_dir, f"chembl_activity_split_{name}.json")
    
    # We pass multiple separate --filter arguments
    cmd = [
        "uv", "run", chembl_script, "activity",
        "--filter", f"target_chembl_id={tid}",
        "--filter", "standard_type=IC50",
        "--limit", "30",
        "--normalize",
        "--output", out_file
    ]
    
    print(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error for {name}: {res.stderr}")
    else:
        # Verify if we got results
        if os.path.exists(out_file):
            with open(out_file, 'r') as f:
                try:
                    data = json.load(f)
                    count = data.get("page_meta", {}).get("total_count", 0)
                    print(f"  Success for {name}: Found {count} records.")
                except Exception as e:
                    print(f"  Error reading output for {name}: {e}")
