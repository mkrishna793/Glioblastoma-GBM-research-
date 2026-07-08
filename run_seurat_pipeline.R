library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(future)

# Enable parallel processing
print("Configuring multicore parallelization (32 workers)...")
plan("multicore", workers = 32)
options(future.globals.maxSize = 24 * 1024^3) # 24 GB limit for globals

print("Loading scRNAseq.rdata...")
load("/home/zeus/scRNAseq.rdata")

print("Initializing Seurat objects...")
seurat_list <- list()
for (i in 1:16) {
  var_name <- paste0("sc", i)
  mat <- get(var_name)
  seurat_list[[i]] <- CreateSeuratObject(counts = mat, project = paste0("Patient", i))
}

print("Merging samples into a single object...")
merged_seurat <- merge(seurat_list[[1]], y = unlist(seurat_list[2:16]), add.cell.ids = paste0("P", 1:16))
print(paste("Merged object cells count:", ncol(merged_seurat)))

# Clean up memory
rm(sc1, sc2, sc3, sc4, sc5, sc6, sc7, sc8, sc9, sc10, sc11, sc12, sc13, sc14, sc15, sc16, seurat_list)
gc()

print("Calculating mitochondrial percentage...")
# Human mitochondrial genes start with "MT-"
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

print("Filtering cells (Quality Control)...")
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
print(paste("Cells remaining after QC:", ncol(merged_seurat)))

print("Normalizing data...")
merged_seurat <- NormalizeData(merged_seurat)

print("Finding variable features...")
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

print("Scaling data...")
merged_seurat <- ScaleData(merged_seurat)

print("Running PCA...")
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat), npcs = 30)

print("Finding neighbors...")
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20)

print("Finding clusters (resolution = 0.4)...")
merged_seurat <- FindClusters(merged_seurat, resolution = 0.4)

print("Running UMAP...")
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)

# Save Seurat object
print("Saving merged Seurat object...")
saveRDS(merged_seurat, file = "/home/zeus/merged_seurat.rds")

# Generate UMAP Plot by Cluster
print("Generating UMAP plot by cluster...")
png("/home/zeus/umap_by_cluster.png", width = 800, height = 800, res = 150)
print(DimPlot(merged_seurat, reduction = "umap", label = TRUE, pt.size = 0.1) + NoLegend() + ggtitle("UMAP Clustering (76,000+ Glioma Cells)"))
dev.off()

# Generate UMAP Plot by Patient
print("Generating UMAP plot by patient...")
png("/home/zeus/umap_by_patient.png", width = 1000, height = 800, res = 150)
print(DimPlot(merged_seurat, reduction = "umap", group.by = "orig.ident", pt.size = 0.1) + ggtitle("Single-Cell UMAP by Patient Cohort"))
dev.off()

# Generate DotPlot of Cell Type Markers
print("Generating cell-type marker DotPlot...")
markers <- c("RBFOX3", "SNAP25", "GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CSPG4", "CD68", "CD14", "PTPRC", "CLDN5", "FLT1")
png("/home/zeus/dotplot_markers.png", width = 1200, height = 600, res = 150)
print(DotPlot(merged_seurat, features = markers) + RotatedAxis() + ggtitle("Cell-Type Annotation Markers"))
dev.off()

# Generate Violin Plots of target pathways (synaptic, glioma-state, loop genes)
print("Generating Violin Plots of target pathways...")
targets <- c("CD44", "OLIG2", "NLGN3", "GRIA2", "STAT3", "EGFR")
png("/home/zeus/vlnplot_targets.png", width = 1200, height = 800, res = 150)
print(VlnPlot(merged_seurat, features = targets, ncol = 3, pt.size = 0) + ggtitle("Expression of Glioblastoma Pathway Targets"))
dev.off()

print("Pipeline completed successfully!")
