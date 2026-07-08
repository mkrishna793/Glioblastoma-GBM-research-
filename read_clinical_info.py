import pandas as pd

file_path = r"C:\Users\bhanu\.gemini\antigravity\scratch\gbm_analysis\CGGA_IDH_A_Clinical_Information.xlsx"
df = pd.read_excel(file_path)
print("Shape:", df.shape)
print("Columns:", list(df.columns))
print("\nFirst 5 rows:")
print(df.head(5).to_string())
print("\nUnique values in key columns:")
for col in ['Histology', 'Grade', 'Gender', 'IDH_mutation_status', '1p19q_codeletion_status', 'OS', 'Censor']:
    if col in df.columns:
        print(f"{col}:", df[col].unique())
