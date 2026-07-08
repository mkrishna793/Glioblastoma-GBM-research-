import pandas as pd
import os

print("Loading matched ClinVar pathogenic variants...")
df_matches = pd.read_csv(r"D:\research of the GBM\clinvar_pathogenic_matches.txt", sep="\t")

print("Loading GLASS recovered mutations...")
df_barcodes = pd.read_csv(r"D:\new data of the GBM research\GLASS_recovered_mutations.txt", sep="\t")

# Merge to get patient barcodes
df_patient_matches = pd.merge(df_matches, df_barcodes[['variant_id', 'case_barcode']], on='variant_id', how='inner')
print(f"Total patient-mutation entries matching ClinVar pathogenic: {len(df_patient_matches)}")

# Map patient barcodes to RNA-seq patient IDs
# In variants, case_barcode looks like "GLSS-19-0271"
# In RNA-seq, the column looks like "GLSS-19-0271-TP" or "GLSS-19-0271-R1" or similar
# Let's extract the patient ID as the first three parts (e.g., GLSS-19-0271)
df_patient_matches['patient_id'] = df_patient_matches['case_barcode']

# Load RNA-seq sample names
df_rna = pd.read_csv(r"D:\new data of the GBM research\gene_tpm_matrix_all_samples.tsv", sep="\t", nrows=5)
rna_samples = df_rna.columns[1:].tolist() # skip index column

rna_patients = []
for col in rna_samples:
    parts = col.replace("-", ".").split(".")
    pid = "-".join(parts[:3]) # Standardize to GLSS-XX-XXXX
    rna_patients.append(pid)

rna_patients = set(rna_patients)
print(f"Total unique patients in RNA-seq: {len(rna_patients)}")

# Check overlap
variant_patients = set(df_patient_matches['patient_id'].unique())
overlap_patients = rna_patients.intersection(variant_patients)
print(f"Patients overlapping between RNA-seq and ClinVar pathogenic mutations: {len(overlap_patients)}")

# Breakdown of overlap by gene
df_overlap = df_patient_matches[df_patient_matches['patient_id'].isin(overlap_patients)]
print("\nBreakdown of patients carrying ClinVar pathogenic mutations in the RNA-seq cohort:")
for gene in ['IDH1', 'ATRX', 'TP53']:
    gene_patients = df_overlap[df_overlap['gene_symbol'] == gene]['patient_id'].unique()
    print(f"  {gene}: {len(gene_patients)} patients")
    print(f"    Patient IDs: {gene_patients}")
