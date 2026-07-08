import pandas as pd
import os

data_dir = r"D:\new data of the GBM research"
out_path = r"D:\research of the GBM\true_mutation_counts.txt"

print("Auditing unique mutations per aliquot...")
aliquot_mutation_sets = {}

try:
    reader = pd.read_csv(os.path.join(data_dir, "variants.passgeno.csv.gz"), chunksize=200000)
    for chunk in reader:
        # Group by aliquot and get unique coordinates (chrom, pos, alt)
        for name, group in chunk.groupby('aliquot_barcode'):
            # Convert to a set of tuples
            coords = set(zip(group['chrom'], group['pos'], group['alt']))
            if name not in aliquot_mutation_sets:
                aliquot_mutation_sets[name] = set()
            aliquot_mutation_sets[name].update(coords)
except EOFError:
    pass

# Compile counts
counts_data = []
for name, coords in aliquot_mutation_sets.items():
    counts_data.append({
        'aliquot_barcode': name,
        'true_mutation_count': len(coords)
    })

df_true_counts = pd.DataFrame(counts_data)
df_true_counts['case_barcode'] = df_true_counts['aliquot_barcode'].apply(lambda x: "-".join(x.split("-")[:3]) if isinstance(x, str) else "")
df_true_counts['timepoint'] = df_true_counts['aliquot_barcode'].apply(lambda x: x.split("-")[3] if isinstance(x, str) and len(x.split("-")) > 3 else "UNK")

print(f"Total aliquots processed: {len(df_true_counts)}")
print("\nTop 10 aliquots by true mutation count:")
print(df_true_counts.sort_values(by='true_mutation_count', ascending=False).head(10).to_string())

# Save to file
df_true_counts.to_csv(out_path, sep="\t", index=False)
print(f"True mutation counts saved to {out_path}")
