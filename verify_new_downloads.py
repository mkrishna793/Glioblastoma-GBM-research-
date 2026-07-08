import gzip
import os

data_dir = r"D:\new data of the GBM research"
gz_files = ['variants.anno.csv.gz', 'variants.passgeno.csv.gz']

print("Verifying integrity of newly downloaded files...")
for f in gz_files:
    path = os.path.join(data_dir, f)
    if not os.path.exists(path):
        print(f"  {f}: NOT FOUND")
        continue
    size_mb = os.path.getsize(path) / (1024*1024)
    print(f"  {f}: Size is {size_mb:.2f} MB. Reading file to the end to check for corruption...")
    try:
        count = 0
        with gzip.open(path, 'rt') as fh:
            for line in fh:
                count += 1
                if count % 2000000 == 0:
                    print(f"    Read {count} lines...")
        print(f"  {f}: OK. Completed reading all {count} lines successfully. No corruption found.")
    except Exception as e:
        print(f"  {f}: CORRUPTED/INCOMPLETE - {type(e).__name__}: {e}")
