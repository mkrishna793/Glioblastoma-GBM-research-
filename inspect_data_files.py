import zipfile
import os
import io
import pandas as pd

base_dir = r"D:\Glioblastoma Multiforme GBM new research files"
zip_files = [
    "CGGA_IDH_A_Methylation_EPIC_Array_20250915.zip",
    "CGGA_IDH_A_RNAseq_RSEM_20250915.zip",
    "CGGA_IDH_A_Phosphoproteomics_MS_Abundance_20250915.zip",
    "CGGA_IDH_A_Proteomics_MS_Abundance_20250915.zip"
]

for z_name in zip_files:
    z_path = os.path.join(base_dir, z_name)
    print("========================================")
    print("Zip File:", z_name)
    print("========================================")
    with zipfile.ZipFile(z_path, 'r') as z:
        for name in z.namelist():
            if name.endswith('.txt'):
                print(f"Reading first few lines of {name}...")
                with z.open(name) as f:
                    # Read first 3 lines
                    wrapper = io.TextIOWrapper(f, encoding='utf-8')
                    for i in range(3):
                        line = wrapper.readline()
                        if not line:
                            break
                        # Print first 200 characters of each line to avoid clutter
                        print(f"Line {i+1}: {line[:200]}...")
                
                # Check shape by reading the file headers and counting lines/columns
                # Since we don't want to load 435MB of methylation in memory, let's read it chunk-wise or just show column count
                with z.open(name) as f:
                    wrapper = io.TextIOWrapper(f, encoding='utf-8')
                    header = wrapper.readline()
                    sample_ids = header.strip().split('\t')
                    print(f"Number of columns (samples/features): {len(sample_ids)}")
                    print(f"First 10 column names: {sample_ids[:10]}")
