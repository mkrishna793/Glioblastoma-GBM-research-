import pandas as pd
import os

data_dir = r"D:\new data of the GBM research"
target_aliquots = [
    'GLSS-HK-0005-TP-01D-WXS-0GCSJ3.0',
    'GLSS-HK-0005-TP-01D-WXS-0GCSJ3.1',
    'TCGA-14-1402-TP-01D-WXS-4A9FK1.0',
    'GLSS-HK-0005-R1-01D-WXS-123456' # let's see any recurrence
]

print("Auditing passgeno file for target aliquots...")
extracted_rows = []

try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=200000)
    for chunk in reader:
        # Check if any target aliquots are in this chunk
        sub = chunk[chunk['aliquot_barcode'].isin(target_aliquots)]
        if not sub.empty:
            extracted_rows.append(sub)
except EOFError:
    pass

if extracted_rows:
    df_sub = pd.concat(extracted_rows, ignore_index=True)
    print(f"Total rows extracted for target aliquots: {len(df_sub)}")
    print("\nSample rows:")
    print(df_sub.head(15).to_string())
    print("\nValue counts of aliquot_barcodes:")
    print(df_sub['aliquot_barcode'].value_counts())
    print("\nUnique variants in each aliquot:")
    for al in df_sub['aliquot_barcode'].unique():
        print(f"  {al}: {len(df_sub[df_sub['aliquot_barcode'] == al]['variant_id'].unique())} unique variants")
else:
    print("No rows found for target aliquots.")
