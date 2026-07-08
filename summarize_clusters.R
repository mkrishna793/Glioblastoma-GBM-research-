library(Seurat)
library(Matrix)

# Load the merged Seurat object
merged_seurat <- readRDS("/home/zeus/merged_seurat.rds")

# Print target genes average expression
targets <- c("CD44", "OLIG2", "NLGN3", "GRIA2", "STAT3", "EGFR")
avg_exp <- AverageExpression(merged_seurat, features = targets, slot = "data")$RNA
avg_exp_t <- t(avg_exp)

print("--- TARGET GENES MEAN EXPRESSION PER CLUSTER ---")
write.csv(avg_exp_t, file = "/home/zeus/target_genes_cluster_summary.csv")
print(avg_exp_t)

# Print cell-type markers average expression
markers <- c("RBFOX3", "SNAP25", "GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CSPG4", "CD68", "CD14", "CLDN5", "FLT1")
avg_markers <- AverageExpression(merged_seurat, features = markers, slot = "data")$RNA
avg_markers_t <- t(avg_markers)

print("--- CELL-TYPE MARKERS MEAN EXPRESSION PER CLUSTER ---")
write.csv(avg_markers_t, file = "/home/zeus/cell_type_markers_cluster_summary.csv")
print(avg_markers_t)
