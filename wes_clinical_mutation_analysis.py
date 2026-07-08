import pandas as pd
import numpy as np
import os

# Define paths
data_dir = r"D:\research of the GBM"
clinical_file = os.path.join(data_dir, "CGGA.WEseq_286_clinical.20200506.txt")
mutation_file = os.path.join(data_dir, "CGGA.WEseq_286.20200506.txt")
output_report = os.path.join(data_dir, "wes_mutation_clinical_report.txt")

print("Loading clinical data...")
# Read clinical info and strip whitespaces from column headers
df_clin = pd.read_csv(clinical_file, sep="\t")
df_clin.columns = df_clin.columns.str.strip()
print(f"Clinical shape: {df_clin.shape}")
print("Cleaned Columns:", list(df_clin.columns))

print("Loading mutation data...")
# Read WES mutations (sample IDs as columns, genes as index/rows)
df_mut = pd.read_csv(mutation_file, sep="\t", index_col=0, low_memory=False)
df_mut.columns = df_mut.columns.str.strip()
print(f"Mutation shape: {df_mut.shape}")

# Transpose mutations so patients are rows and genes are columns
df_mut_t = df_mut.transpose()
# Fill NaN/empty strings with 0, and non-empty with 1 (binary mutation matrix)
df_mut_bin = df_mut_t.notna().astype(int)

# Create report buffer
report = []
report.append("=========================================")
report.append("CGGA 286 WES CLINICAL & GENOMIC MUTATION REPORT")
report.append("=========================================\n")

# Cohort Summary
report.append("--- 1. CLINICAL COHORT SUMMARY ---")
report.append(f"Total Patients: {len(df_clin)}")
report.append("\nSubtype Distribution:")
report.append(df_clin['Subtype'].value_counts().to_string())
report.append("\nWHO Grade Distribution:")
report.append(df_clin['Grade'].value_counts().to_string())
report.append("\nGender Distribution:")
report.append(df_clin['Gender'].value_counts().to_string())
report.append("\nChemotherapy Status (TMZ):")
report.append(df_clin['Chemo_status (TMZ treated=1;un-treated=0)'].value_counts().to_string())
report.append("\nRadiotherapy Status:")
report.append(df_clin['Radio_status (treated=1;un-treated=0)'].value_counts().to_string())
report.append("\nIDH Mutation Status:")
report.append(df_clin['IDH_mut_status'].value_counts().to_string())
report.append("\nMGMT Promoter Methylation Status:")
report.append(df_clin['MGMTp_methylation_status'].value_counts().to_string())

# Mutation Frequency Analysis
report.append("\n\n--- 2. TOP SOMATIC MUTATED GENES ---")
mutation_rates = df_mut_bin.mean() * 100
top_mutated = mutation_rates.sort_values(ascending=False).head(30)
report.append("Gene Name\tMutation Rate (%)")
for gene, rate in top_mutated.items():
    report.append(f"{gene}\t{rate:.2f}%")

# Check specific key genes
report.append("\n\n--- 3. MUTATION RATES OF TARGET GBM DRIVERS ---")
target_drivers = ['IDH1', 'TP53', 'ATRX', 'PTEN', 'EGFR', 'CIC', 'OBSCN', 'AHNAK2', 'TERT']
for g in target_drivers:
    if g in df_mut_bin.columns:
        rate = df_mut_bin[g].mean() * 100
        report.append(f"{g}: {rate:.2f}%")
    else:
        report.append(f"{g}: 0.00% (Not found/not mutated)")

# Connect WES Mutations to Clinical Survival Outcomes
# Merge clinical and binary mutation matrices on Patient ID
df_clin.set_index('CGGA_ID', inplace=True)
# Ensure indexes are string and stripped
df_clin.index = df_clin.index.astype(str).str.strip()
df_mut_bin.index = df_mut_bin.index.astype(str).str.strip()

df_merged = df_clin.join(df_mut_bin, how='inner')
print(f"Merged WES shape: {df_merged.shape}")

# Convert OS and Censor columns to numeric, drop missing values
df_merged['OS'] = pd.to_numeric(df_merged['OS'], errors='coerce')
df_merged['Censor (alive=0; dead=1)'] = pd.to_numeric(df_merged['Censor (alive=0; dead=1)'], errors='coerce')
df_valid_os = df_merged.dropna(subset=['OS', 'Censor (alive=0; dead=1)'])

report.append("\n\n--- 4. SOMATIC MUTATIONS & PATIENT SURVIVAL CORRELATION ---")

# Let's compare OBSCN mutation vs. Wildtype survival
for g in ['OBSCN', 'AHNAK2', 'TP53', 'PTEN', 'IDH1']:
    if g in df_valid_os.columns:
        mut_group = df_valid_os[df_valid_os[g] == 1]
        wt_group = df_valid_os[df_valid_os[g] == 0]
        
        mut_mean_os = mut_group['OS'].mean()
        wt_mean_os = wt_group['OS'].mean()
        
        report.append(f"\n{g} Mutation Status:")
        report.append(f"  Mutated Patients count: {len(mut_group)}")
        report.append(f"  Mutated Group Mean OS: {mut_mean_os:.1f} days")
        report.append(f"  Wildtype Group Mean OS: {wt_mean_os:.1f} days")
        
        # Calculate Censor details
        mut_death_rate = mut_group['Censor (alive=0; dead=1)'].mean() * 100
        wt_death_rate = wt_group['Censor (alive=0; dead=1)'].mean() * 100
        report.append(f"  Mutated Group Death Rate: {mut_death_rate:.1f}%")
        report.append(f"  Wildtype Group Death Rate: {wt_death_rate:.1f}%")

# Save report
report_text = "\n".join(report)
with open(output_report, 'w') as f:
    f.write(report_text)

print("Report successfully written to:", output_report)
