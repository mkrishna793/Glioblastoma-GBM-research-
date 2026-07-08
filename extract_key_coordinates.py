import pandas as pd
import os

data_dir = r"D:\new data of the GBM research"
out = r"C:\Users\bhanu\.gemini\antigravity\brain\774b6382-5a3e-4f40-985d-ba97c732d8ea\key_genes_coordinates.txt"

print("Extracting coordinates for key genes (IDH1, ATRX, TP53)...")
key_genes = ['IDH1', 'ATRX', 'TP53']
extracted = []

try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.anno.csv.gz"), chunksize=200000)
    for chunk in reader:
        # Filter for key genes and coding mutations
        subset = chunk[
            chunk['gene_symbol'].isin(key_genes) & 
            chunk['variant_classification'].isin([
                'MISSENSE', 'NONSENSE', 'FRAME_SHIFT_DEL', 'FRAME_SHIFT_INS',
                'SPLICE_SITE', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'NONSTOP'
            ])
        ][['variant_id', 'gene_symbol', 'chrom', 'pos', 'ref', 'alt', 'variant_classification']]
        extracted.append(subset)
except EOFError:
    print("  Reached EOF of variants.anno file.")
except Exception as e:
    print(f"  Error: {e}")

if extracted:
    df_extracted = pd.concat(extracted, ignore_index=True)
    print(f"Total mutations extracted: {len(df_extracted)}")
    df_extracted.to_csv(out, sep="\t", index=False)
    print(f"Saved to {out}")
else:
    print("No mutations extracted.")
