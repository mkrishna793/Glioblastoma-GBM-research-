import pandas as pd
import gzip
import os

data_dir = r"D:\new data of the GBM research"
out = r"D:\new data of the GBM research\GLASS_recovered_mutations.txt"

print("Salvaging coding mutations from truncated variants.anno.csv.gz using pandas chunks...")
chunks = []
try:
    # Use pandas read_csv in chunks
    reader = pd.read_csv(os.path.join(data_dir, "variants.anno.csv.gz"), chunksize=100000)
    for chunk in reader:
        # Filter for coding mutations in current chunk (all-caps)
        coding_chunk = chunk[chunk['variant_classification'].isin([
            'MISSENSE', 'NONSENSE', 'FRAME_SHIFT_DEL', 'FRAME_SHIFT_INS',
            'SPLICE_SITE', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'NONSTOP'
        ])][['variant_id', 'gene_symbol', 'variant_classification']]
        chunks.append(coding_chunk)
except EOFError:
    print("  Reached EOF marker in truncated variants.anno file.")
except Exception as e:
    print(f"  Error reading variants.anno: {e}")

if chunks:
    df_coding = pd.concat(chunks, ignore_index=True)
    print(f"Salvaged coding mutations: {len(df_coding)}")
else:
    df_coding = pd.DataFrame()
    print("No coding mutations salvaged.")

print("\nSalvaging pass genotypes from truncated variants.passgeno.csv.gz using pandas chunks...")
geno_chunks = []
try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=200000)
    for chunk in reader:
        # Filter for pass call = 1
        pass_chunk = chunk[chunk['ssm2_pass_call'] == 1][['case_barcode', 'variant_id']]
        geno_chunks.append(pass_chunk)
except EOFError:
    print("  Reached EOF marker in truncated genotypes file.")
except Exception as e:
    print(f"  Error reading genotypes: {e}")

if geno_chunks:
    df_geno = pd.concat(geno_chunks, ignore_index=True)
    print(f"Salvaged genotypes: {len(df_geno)}")
else:
    df_geno = pd.DataFrame()
    print("No genotypes salvaged.")

# Merge and print counts
if not df_coding.empty and not df_geno.empty:
    df_merged = pd.merge(df_geno, df_coding, on='variant_id')
    key_genes = ['IDH1', 'ATRX', 'TP53', 'PTEN', 'EGFR']
    print("\nSalvaged Mutation Counts:")
    for g in key_genes:
        pts = df_merged[df_merged['gene_symbol'] == g]['case_barcode'].unique()
        print(f"  {g}-mutant patients: {len(pts)}")
    
    df_merged.to_csv(out, sep="\t", index=False)
    print(f"\nSaved salvaged mutations to {out}")
else:
    print("\nFailed to salvage enough data.")
