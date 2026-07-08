load("/home/zeus/scRNAseq.rdata")
objs <- ls()
print("Loaded objects:")
print(objs)
for (obj in objs) {
  print(paste("Class of", obj, ":"))
  print(class(get(obj)))
}
