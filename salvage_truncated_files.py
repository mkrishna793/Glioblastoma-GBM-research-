import pandas as pd
import numpy as np
import gzip
import os, sys

sys.stdout.reconfigure(encoding='utf-8')
data_dir = r"D:\new data of the GBM research"
out = r"D:\new data of the GBM research\GLASS_recovered_mutations.txt"

print("Salvaging coding mutations from truncated variants.anno.csv.gz...")
coding_variants = []
columns = []

try:
    with gzip.open(os.path.join(data_dir, "variants.anno.csv.gz"), "rt") as f:
        # Read header
        header = f.readline().strip()
        columns = header.split(",")
        
        # Read lines one by one, catching EOFError
        count = 0
        while True:
            try:
                line = f.readline()
                if not line:
                    break
                
                parts = line.strip().split(",")
                if len(parts) >= len(columns):
                    # We look for coding mutations in variant_classification (index 6)
                    classification = parts[6]
                    if classification in ['Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'In_Frame_Del', 'In_Frame_Ins', 'Nonstop_Mutation']:
                        # Save mapping of variant_id to gene_symbol
                        coding_variants.append({
                            'variant_id': int(parts[0]),
                            'gene_symbol': parts[5],
                            'variant_classification': classification
                        })
                count += 1
                if count % 500000 == 0:
                    print(f"  Processed {count} lines, salvaged {len(coding_variants)} coding mutations...")
            except EOFError:
                print("  Reached EOF marker in truncated file. Stop reading.")
                break
except Exception as e:
    print(f"  Error reading file: {e}")

print(f"Completed salvaging. Total coding mutations recovered: {len(coding_variants)}")
df_coding = pd.DataFrame(coding_variants)

# Now salvage genotypes
print("\nSalvaging genotypes from truncated variants.passgeno.csv.gz...")
genotypes = []
try:
    with gzip.open(os.path.join(data_dir, "variants.passgeno.csv.gz"), "rt") as f:
        header = f.readline().strip()
        geno_cols = header.split(",")
        
        count = 0
        while True:
            try:
                line = f.readline()
                if not line:
                    break
                
                parts = line.strip().split(",")
                if len(parts) >= len(geno_cols):
                    # Index 9 is ssm2_pass_call (usually 1.0 or 1)
                    pass_call = parts[9]
                    if pass_call in ['1', '1.0', 1, 1.0]:
                        genotypes.append({
                            'case_barcode': parts[2],
                            'variant_id': int(parts[1])
                        })
                count += 1
                if count % 500000 == 0:
                    print(f"  Processed {count} lines, salvaged {len(genotypes)} pass genotypes...")
            except EOFError:
                print("  Reached EOF marker in truncated genotypes file. Stop reading.")
                break
except Exception as e:
    print(f"  Error reading genotypes: {e}")

print(f"Completed salvaging genotypes. Total pass genotypes recovered: {len(genotypes)}")
df_geno = pd.DataFrame(genotypes)

# Match mutations to patients
if not df_coding.empty and not df_geno.empty:
    df_merged = pd.merge(df_geno, df_coding, on='variant_id')
    
    # Let's count unique patients with key mutations
    key_genes = ['IDH1', 'ATRX', 'TP53', 'PTEN', 'EGFR']
    for g in key_genes:
        pts = df_merged[df_merged['gene_symbol'] == g]['case_barcode'].unique()
        print(f"  {g}-mutant patients: {len(pts)}")
        
    # Save the salvaged mutation mapping
    df_merged.to_csv(out, sep="\t", index=False)
    print(f"Salvaged mutation map saved to {out}")
else:
    print("Failed to salvage enough records.")
