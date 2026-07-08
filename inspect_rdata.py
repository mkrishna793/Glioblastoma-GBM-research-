import pyreadr
try:
    result = pyreadr.read_r(r"D:\Glioblastoma Multiforme GBM new research files\scRNAseq.rdata")
    print("Objects found in RData:", list(result.keys()))
except Exception as e:
    print("Error reading RData file using pyreadr:", e)
