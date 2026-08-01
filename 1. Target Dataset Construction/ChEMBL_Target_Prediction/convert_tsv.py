import pandas as pd

# 1. Load and Clean ChEMBL Targets
# Your file actually contains UniProt IDs in the 'Accessions' column!
print("Reading ChEMBL Targets...")
chembl_df = pd.read_csv("targets_tsv.tsv", sep='\t')

# Extract UniProt IDs (handling multiple IDs like 'P04629|Q16620')
chembl_uniprots = set()
for accs in chembl_df['Accessions'].dropna():
    for acc in accs.split('|'):
        chembl_uniprots.add(acc.strip())

print(f"✅ Extracted {len(chembl_uniprots)} unique UniProt IDs from ChEMBL.")

# 2. Load Swiss Targets
swiss_df = pd.read_csv("SwissTargetPrediction.txt")
swiss_uniprots = set(swiss_df['Uniprot ID'].unique())
print(f"✅ Extracted {len(swiss_uniprots)} unique UniProt IDs from Swiss.")

# 3. Create 'Predicted Pool' (Union of AI tools)
predicted_targets = chembl_uniprots.union(swiss_uniprots)
print(f"Total Unique Predicted Targets: {len(predicted_targets)}")

# 4. Prepare HL-60 Gene List for Mapping
# We need to convert these Gene Symbols (NRAS, etc.) to UniProt to do the intersection
hl60_df = pd.read_csv("HL_60_RNAi.txt", header=None, names=['Gene'])
hl60_genes = hl60_df['Gene'].unique()

print(f"\n✅ Found {len(hl60_genes)} Essential Genes in HL-60.")
print("⬇️ SAVE THIS LIST FOR UNIPROT MAPPING ⬇️")
hl60_df.to_csv("HL60_Genes_for_Mapping.txt", index=False, header=False)
print("File 'HL60_Genes_for_Mapping.txt' saved.")