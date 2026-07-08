load("/home/zeus/scRNAseq.rdata")
for (i in 1:16) {
  mat <- get(paste0("sc", i))
  print(paste("Sample", i, "dimensions:", paste(dim(mat), collapse=" x ")))
}
