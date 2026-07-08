import pandas as pd
import os

data_dir = r"D:\new data of the GBM research"
out = r"D:\research of the GBM\key_variants_with_barcodes.txt"

print("Loading matched ClinVar pathogenic variants...")
df_matches = pd.read_csv(r"D:\research of the GBM\clinvar_pathogenic_matches.txt", sep="\t")
target_ids = set(df_matches['variant_id'].unique())
print(f"Targeting {len(target_ids)} variant IDs...")

# Read variants.passgeno.csv.gz in chunks and match variant_id
matched_barcodes = []
try:
    # Let's inspect the columns first to see if it is aliquot_barcode or case_barcode
    first_chunk = next(pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=5))
    cols = first_chunk.columns.tolist()
    print("Columns in passgeno:", cols)
    
    barcode_col = 'aliquot_barcode' if 'aliquot_barcode' in cols else ('case_barcode' if 'case_barcode' in cols else cols[0])
    
    reader = pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=200000)
    for chunk in reader:
        subset = chunk[chunk['variant_id'].isin(target_ids)][['variant_id', barcode_col, 'af']]
        matched_barcodes.append(subset)
except EOFError:
    print("  Reached EOF of passgeno file.")
except Exception as e:
    print(f"  Error: {e}")

if matched_barcodes:
    df_barcodes = pd.concat(matched_barcodes, ignore_index=True)
    # Standardize barcode_col to case_barcode (first 3 parts)
    # aliquot_barcode looks like "GLSS-19-0268-TP-01D-WXS-12345"
    # case_barcode looks like "GLSS-19-0268"
    df_barcodes['case_barcode'] = df_barcodes[barcode_col].apply(lambda x: "-".join(x.split("-")[:3]) if isinstance(x, str) else "")
    print(f"Matched {len(df_barcodes)} variant-aliquot links.")
    
    print("df_barcodes columns before merge:", df_barcodes.columns.tolist())
    # Merge back with matches
    df_final = pd.merge(df_matches, df_barcodes, on='variant_id', how='inner')
    print("df_final columns:", df_final.columns.tolist())
    print(f"Merged successfully: {len(df_final)} rows.")
    df_final.to_csv(out, sep="\t", index=False)
    print(f"Saved to {out}")
else:
    print("No matches found in passgeno.")
