import pandas as pd
import numpy as np
import os

data_dir = r"D:\new data of the GBM research"
out_path = r"D:\research of the GBM\hypothesis_4_results.txt"

print("Extracting MMR variants from variants.anno...")
mmr_genes = ['MSH6', 'MSH2', 'MLH1', 'PMS2']
mmr_extracted = []

try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.anno.csv.gz"), chunksize=200000)
    for chunk in reader:
        subset = chunk[
            chunk['gene_symbol'].isin(mmr_genes) & 
            chunk['variant_classification'].isin([
                'MISSENSE', 'NONSENSE', 'FRAME_SHIFT_DEL', 'FRAME_SHIFT_INS',
                'SPLICE_SITE', 'IN_FRAME_DEL', 'IN_FRAME_INS', 'NONSTOP'
            ])
        ][['variant_id', 'gene_symbol', 'chrom', 'pos', 'ref', 'alt', 'variant_classification']]
        mmr_extracted.append(subset)
except EOFError:
    pass
except Exception as e:
    print(f"  Error: {e}")

if mmr_extracted:
    df_mmr_vars = pd.concat(mmr_extracted, ignore_index=True)
    print(f"Total coding MMR variants extracted: {len(df_mmr_vars)}")
else:
    df_mmr_vars = pd.DataFrame()
    print("No MMR variants extracted.")

# Get patient barcodes for these MMR variants from passgeno
df_mmr_patients = pd.DataFrame()
if not df_mmr_vars.empty:
    mmr_ids = set(df_mmr_vars['variant_id'].unique())
    matched_mmr = []
    try:
        reader = pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=200000)
        for chunk in reader:
            subset = chunk[chunk['variant_id'].isin(mmr_ids)][['variant_id', 'aliquot_barcode', 'af']]
            matched_mmr.append(subset)
    except EOFError:
        pass
        
    if matched_mmr:
        df_mmr_barcodes = pd.concat(matched_mmr, ignore_index=True)
        df_mmr_barcodes['case_barcode'] = df_mmr_barcodes['aliquot_barcode'].apply(lambda x: "-".join(x.split("-")[:3]) if isinstance(x, str) else "")
        df_mmr_barcodes['timepoint'] = df_mmr_barcodes['aliquot_barcode'].apply(lambda x: x.split("-")[3] if isinstance(x, str) and len(x.split("-")) > 3 else "UNK")
        
        df_mmr_patients = pd.merge(df_mmr_vars, df_mmr_barcodes, on='variant_id', how='inner')
        print(f"Matched {len(df_mmr_patients)} patient-MMR variant entries.")

# Calculate total somatic mutation count per aliquot from passgeno
print("Calculating mutation counts per aliquot...")
mutation_counts = {}
try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=200000)
    for chunk in reader:
        counts = chunk['aliquot_barcode'].value_counts()
        for barcode, count in counts.items():
            mutation_counts[barcode] = mutation_counts.get(barcode, 0) + count
except EOFError:
    pass

df_counts = pd.DataFrame(list(mutation_counts.items()), columns=['aliquot_barcode', 'mutation_count'])
df_counts['case_barcode'] = df_counts['aliquot_barcode'].apply(lambda x: "-".join(x.split("-")[:3]) if isinstance(x, str) else "")
df_counts['timepoint'] = df_counts['aliquot_barcode'].apply(lambda x: x.split("-")[3] if isinstance(x, str) and len(x.split("-")) > 3 else "UNK")

print(f"Total aliquots counted: {len(df_counts)}")

# Write results
with open(out_path, "w") as f:
    f.write("="*80 + "\n")
    f.write("HYPOTHESIS 4 TESTING: DNA DAMAGE CHECKPOINT ESCAPE VIA MMR HYPERMUTATION\n")
    f.write("="*80 + "\n\n")
    
    if not df_mmr_patients.empty:
        f.write("1. Coding Mismatch Repair (MMR) variants in the cohort:\n")
        f.write(df_mmr_patients[['case_barcode', 'timepoint', 'gene_symbol', 'variant_classification', 'af']].to_string() + "\n\n")
        
        # Check for hypermutation (defined as mutation_count > 1000 in WES)
        f.write("2. Mutation counts of patients with MMR mutations:\n")
        df_mmr_counts = pd.merge(df_mmr_patients, df_counts, on=['aliquot_barcode', 'case_barcode', 'timepoint'], how='inner')
        f.write(df_mmr_counts[['case_barcode', 'timepoint', 'gene_symbol', 'variant_classification', 'mutation_count']].to_string() + "\n\n")
        
        # Check overall hypermutation rate in cohort
        hypermutated = df_counts[df_counts['mutation_count'] >= 1000]
        f.write(f"3. Hypermutation burden in the cohort:\n")
        f.write(f"  Aliquots with >= 1000 mutations (hypermutated): {len(hypermutated)}\n")
        f.write(hypermutated[['aliquot_barcode', 'case_barcode', 'timepoint', 'mutation_count']].to_string() + "\n\n")
        
        # Overlap between hypermutation and MMR mutations
        mmr_aliquots = set(df_mmr_patients['aliquot_barcode'].unique())
        hypermut_aliquots = set(hypermutated['aliquot_barcode'].unique())
        overlap = mmr_aliquots.intersection(hypermut_aliquots)
        f.write(f"4. Overlap between MMR mutations and hypermutation:\n")
        f.write(f"  MMR mutant aliquots that are hypermutated: {len(overlap)} / {len(mmr_aliquots)}\n")
        if overlap:
            f.write(f"    Aliquots: {overlap}\n")
    else:
        f.write("No MMR variants found.\n")

print(f"Results saved to {out_path}")
# Print the contents of the file directly to stdout for verification
with open(out_path, "r") as f:
    print("\n--- FILE CONTENTS ---")
    print(f.read())
    print("--- END OF FILE ---")
if not df_mmr_patients.empty:
    print(f"Unique MMR mutant patients: {len(df_mmr_patients['case_barcode'].unique())}")
    hyper = df_counts[df_counts['mutation_count'] >= 1000]
    print(f"Hypermutated aliquots in cohort: {len(hyper)}")
