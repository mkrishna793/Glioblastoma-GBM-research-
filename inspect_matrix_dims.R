load("/home/zeus/scRNAseq.rdata")
objs <- paste0("sc", 1:16)
total_cells <- 0
for (obj in objs) {
  mat <- get(obj)
  dims <- dim(mat)
  print(paste("Sample:", obj, "| Genes:", dims[1], "| Cells:", dims[2]))
  total_cells <- total_cells + dims[2]
  
  # Check first 5 genes
  print("First 5 genes:")
  print(head(rownames(mat), 5))
}
print(paste("Total cells across all samples:", total_cells))
