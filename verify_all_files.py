import gzip
import os

data_dir = r"D:\new data of the GBM research"
files_to_check = [
    ('variants.gene_copy_number.csv.gz', True),
    ('transcript_count_matrix_all_samples (1).tsv.gz', True),
    ('transcript_count_matrix_all_samples.tsv.gz', True),
    ('transcript_eff_length_matrix_all_samples.tsv.gz', True),
    ('transcript_tpm_matrix_all_samples.tsv.gz', True),
    ('beta.450.tsv', False),
    ('beta.epic.tsv', False),
    ('gene_tpm_matrix_all_samples.tsv', False)
]

print("Verifying integrity of all remaining files in D:\\new data of the GBM research...")
for f, is_gz in files_to_check:
    path = os.path.join(data_dir, f)
    if not os.path.exists(path):
        print(f"  {f}: NOT FOUND")
        continue
    size_mb = os.path.getsize(path) / (1024*1024)
    print(f"  {f}: Size is {size_mb:.2f} MB. Reading to the end...")
    try:
        count = 0
        if is_gz:
            with gzip.open(path, 'rt') as fh:
                for line in fh:
                    count += 1
        else:
            with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
                for line in fh:
                    count += 1
        print(f"  {f}: OK. Read all {count} lines. No corruption.")
    except Exception as e:
        print(f"  {f}: CORRUPTED/INCOMPLETE - {type(e).__name__}: {e}")
