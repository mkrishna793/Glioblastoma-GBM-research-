import zipfile
import os
import io

base_dir = r"D:\Glioblastoma Multiforme GBM new research files"
zip_files = [
    "CGGA.normal_20.Read_Counts-genes.20230104.txt.zip",
    "CGGA.WEseq_286.20200506.txt.zip",
    "CGGA.WEseq_286_clinical.20200506.txt.zip"
]

for z_name in zip_files:
    z_path = os.path.join(base_dir, z_name)
    print("========================================")
    print("Zip File:", z_name)
    print("========================================")
    with zipfile.ZipFile(z_path, 'r') as z:
        for name in z.namelist():
            print(f"File inside zip: {name}")
            if name.endswith('.txt') or name.endswith('.csv') or name.endswith('.xlsx'):
                print(f"Reading first few lines of {name}...")
                with z.open(name) as f:
                    wrapper = io.TextIOWrapper(f, encoding='utf-8', errors='ignore')
                    for i in range(5):
                        line = wrapper.readline()
                        if not line:
                            break
                        print(f"Line {i+1}: {line[:200].strip()}")
