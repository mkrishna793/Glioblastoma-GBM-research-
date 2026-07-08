import pandas as pd
import numpy as np
import os

# Define paths
data_dir = r"D:\research of the GBM"
normal_file = os.path.join(data_dir, "CGGA.normal_20.Read_Counts-genes.20230104.txt")
tumor_file = os.path.join(data_dir, "CGGA_IDH_A_RNAseq_RSEM_20250915.txt")
output_report = os.path.join(data_dir, "tumor_normal_rna_report.txt")

print("Loading normal brain RNA-Seq data...")
df_normal = pd.read_csv(normal_file, sep="\t", index_col=0)
print(f"Normal RNA shape: {df_normal.shape}")

print("Loading tumor brain RNA-Seq data...")
df_tumor = pd.read_csv(tumor_file, sep="\t", index_col=0)
print(f"Tumor RNA shape: {df_tumor.shape}")

# Normalize normal counts to CPM (Counts Per Million)
print("Normalizing normal counts to CPM...")
total_counts = df_normal.sum(axis=0)
df_normal_cpm = (df_normal / total_counts) * 1e6

# Align genes
common_genes = df_normal_cpm.index.intersection(df_tumor.index)
print(f"Number of common genes: {len(common_genes)}")

# Filter datasets for common genes
df_normal_aligned = df_normal_cpm.loc[common_genes]
df_tumor_aligned = df_tumor.loc[common_genes]

# Calculate means
mean_normal = df_normal_aligned.mean(axis=1)
mean_tumor = df_tumor_aligned.mean(axis=1)

# Calculate Log2 Fold Change (adding a small pseudocount to prevent divide by zero)
log2_fc = np.log2((mean_tumor + 1.0) / (mean_normal + 1.0))

# Create results DataFrame
df_results = pd.DataFrame({
    'Mean_Normal_CPM': mean_normal,
    'Mean_Tumor_RSEM': mean_tumor,
    'Log2_Fold_Change': log2_fc
})

# Find top upregulated genes
upregulated = df_results.sort_values(by='Log2_Fold_Change', ascending=False)

# Build report
report = []
report.append("=========================================")
report.append("CGGA GLIOMA TUMOR VS. NORMAL BRAIN CONTROL RNA-SEQ REPORT")
report.append("=========================================\n")

report.append("--- 1. COMPARISON METRICS ---")
report.append(f"Tumor Samples: {df_tumor.shape[1]}")
report.append(f"Normal Samples: {df_normal.shape[1]}")
report.append(f"Genes Analyzed: {len(common_genes)}")

report.append("\n\n--- 2. TOP 30 UPREGULATED GENES IN TUMOR ---")
report.append("Gene\tMean_Normal\tMean_Tumor\tLog2FC")
for gene, row in upregulated.head(30).iterrows():
    report.append(f"{gene}\t{row['Mean_Normal_CPM']:.2f}\t{row['Mean_Tumor_RSEM']:.2f}\t{row['Log2_Fold_Change']:.2f}")

report.append("\n\n--- 3. TOP 30 DOWNREGULATED GENES IN TUMOR ---")
report.append("Gene\tMean_Normal\tMean_Tumor\tLog2FC")
downregulated = df_results.sort_values(by='Log2_Fold_Change', ascending=True)
for gene, row in downregulated.head(30).iterrows():
    report.append(f"{gene}\t{row['Mean_Normal_CPM']:.2f}\t{row['Mean_Tumor_RSEM']:.2f}\t{row['Log2_Fold_Change']:.2f}")

# Target genes profiles (Subtype markers, synaptic integration, and loops)
target_genes = [
    # Proneural / Stemness
    'OLIG2', 'SOX2', 'ASCL1', 
    # Mesenchymal
    'CD44', 'VIM', 'CHI3L1', 'SPP1',
    # Brain Activity / Synapse integration
    'NLGN3', 'GRIA1', 'GRIA2', 'GRIA3', 'GRIA4',
    # Loops and Drivers
    'EGFR', 'PTEN', 'MGMT', 'STAT3', 'TAZ'
]

report.append("\n\n--- 4. EXPRESSION PROFILE OF TARGET PATHWAY GENES ---")
report.append("Gene\tMean_Normal\tMean_Tumor\tLog2FC\tExpression_Status")
for gene in target_genes:
    # Match case-insensitively or exactly
    found_gene = None
    if gene in df_results.index:
        found_gene = gene
    elif gene.upper() in df_results.index:
        found_gene = gene.upper()
        
    if found_gene:
        row = df_results.loc[found_gene]
        fc = row['Log2_Fold_Change']
        status = "UPREGULATED" if fc > 1.0 else ("DOWNREGULATED" if fc < -1.0 else "NO CHANGE")
        report.append(f"{found_gene}\t{row['Mean_Normal_CPM']:.2f}\t{row['Mean_Tumor_RSEM']:.2f}\t{fc:.2f}\t{status}")
    else:
        report.append(f"{gene}\tN/A\tN/A\tN/A\tNOT FOUND IN COHORT")

# Save report
report_text = "\n".join(report)
with open(output_report, 'w') as f:
    f.write(report_text)

print("Report successfully written to:", output_report)
