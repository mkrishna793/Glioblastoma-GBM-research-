import subprocess
import json
import os

ucsc_script = r"C:\Users\bhanu\.gemini\config\plugins\science\skills\ucsc_conservation_and_tfbs\scripts\get_tfbs.py"
output_dir = r"D:\research of the GBM"

# Coordinates of CHI3L1 promoter region on hg38
# TSS is 203186704. Since negative strand, upstream is 203186704 to 203188704.
coords = "chr1:203186704-203188704"

tracks = ["encRegTfbsClustered", "ReMapTFs"]
results = {}

for track in tracks:
    out_file = os.path.join(output_dir, f"tfbs_{track}_chi3l1.json")
    cmd = [
        "uv", "run", ucsc_script,
        "--coordinates", coords,
        "--tracks", track,
        "--output", out_file
    ]
    
    print(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"  Error for {track}: {res.stderr}")
    else:
        # Parse output file to look for STAT3, NFKB, TEAD, CEBPB
        if os.path.exists(out_file):
            with open(out_file, 'r') as f:
                try:
                    data = json.load(f)
                    # The structure is typically a dictionary containing track results
                    print(f"  Successfully fetched TFBS from {track}")
                    
                    # We will write a small parser inside to extract the hits
                    hits = []
                    # Depending on UCSC JSON structure:
                    # Let's inspect it in the next step or print count
                    items = data.get(track, [])
                    if isinstance(data, dict) and not items:
                        # Sometimes it's nested differently
                        for k, v in data.items():
                            if isinstance(v, list):
                                items = v
                                break
                    
                    print(f"  Found {len(items)} total binding sites in this region.")
                except Exception as e:
                    print(f"  Error reading TFBS file: {e}")
