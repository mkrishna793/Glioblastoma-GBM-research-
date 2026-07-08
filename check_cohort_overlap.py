import pandas as pd
import os

print("Loading matched ClinVar pathogenic variants...")
df_matches = pd.read_csv(r"D:\research of the GBM\clinvar_pathogenic_matches.txt", sep="\t")

print("Loading GLASS recovered mutations (which has patient barcodes)...")
df_barcodes = pd.read_csv(r"D:\new data of the GBM research\GLASS_recovered_mutations.txt", sep="\t")

# Merge matches with barcodes on variant_id
df_patient_matches = pd.merge(df_matches, df_barcodes[['variant_id', 'case_barcode']], on='variant_id', how='inner')
print(f"Total patient-mutation entries matching ClinVar pathogenic: {len(df_patient_matches)}")

print("Loading clinical info for multi-omic cohort...")
clinical_path = r"D:\new data of the GBM research\clinical_info.txt"
if not os.path.exists(clinical_path):
    clinical_path = r"D:\research of the GBM\clinical_info.txt"

if os.path.exists(clinical_path):
    df_clinical = pd.read_csv(clinical_path, sep="\t")
    print(f"Loaded clinical data for {len(df_clinical)} samples.")
    # In clinical_info, the column is case_barcode (looks like GLSS-19-0268, etc)
    clinical_cases = set(df_clinical['case_barcode'].unique())
    variant_cases = set(df_patient_matches['case_barcode'].unique())
    overlap = clinical_cases.intersection(variant_cases)
    print(f"Overlap in case barcodes: {len(overlap)}")
    print(f"Overlapping cases: {overlap}")
    
    # Let's count how many samples of each gene we have in the multi-omic cohort!
    df_overlap_muts = df_patient_matches[df_patient_matches['case_barcode'].isin(overlap)]
    print("\nOverlap mutations details:")
    print(df_overlap_muts[['case_barcode', 'gene_symbol', 'variant_classification', 'clinvar_clinsig']].to_string())
else:
    print("Clinical file not found at expected paths.")
