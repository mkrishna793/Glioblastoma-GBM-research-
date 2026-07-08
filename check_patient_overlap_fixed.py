import pandas as pd

print("Loading final matched variants with barcodes...")
df = pd.read_csv(r"D:\research of the GBM\key_variants_with_barcodes.txt", sep="\t")
print(f"Total entries: {len(df)}")

# Load RNA-seq sample names
df_rna = pd.read_csv(r"D:\new data of the GBM research\gene_tpm_matrix_all_samples.tsv", sep="\t", nrows=5)
rna_samples = df_rna.columns[1:].tolist()

rna_patients = []
for col in rna_samples:
    parts = col.replace("-", ".").split(".")
    pid = "-".join(parts[:3]) # Standardize to GLSS-XX-XXXX
    rna_patients.append(pid)

rna_patients = set(rna_patients)
print(f"Total unique patients in RNA-seq: {len(rna_patients)}")

# Find intersection of patients
variant_patients = set(df['case_barcode'].unique())
overlap_patients = rna_patients.intersection(variant_patients)
print(f"Patients overlapping between RNA-seq and ClinVar pathogenic mutations: {len(overlap_patients)}")

# Breakdown of overlap by gene
df_overlap = df[df['case_barcode'].isin(overlap_patients)]
print("\nBreakdown of patients carrying ClinVar pathogenic mutations in the RNA-seq cohort:")
for gene in ['IDH1', 'ATRX', 'TP53']:
    gene_patients = df_overlap[df_overlap['gene_symbol'] == gene]['case_barcode'].unique()
    print(f"  {gene}: {len(gene_patients)} patients")
    print(f"    Patient IDs: {gene_patients}")
