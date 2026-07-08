import pandas as pd
import os

data_dir = r"D:\research of the GBM"
df_idh_clin = pd.read_excel(os.path.join(data_dir, "CGGA_IDH_A_Clinical_Information.xlsx"))
print("IDH-A columns:", df_idh_clin.columns.tolist())
