import pandas as pd

# Input file from UniProt
input_file = "idmapping_2025_12_10.xlsx - Sheet0.csv"
output_file = "Tier_1_Cleaned_Targets.tsv"

print(f"Reading {input_file}...")
df = pd.read_csv(input_file)

# 1. Create Scoring Columns for Sorting
# Reviewed: 1 if yes, 0 if no
df['is_reviewed'] = df['Reviewed'].apply(lambda x: 1 if x == 'reviewed' else 0)
# Has PDB: 1 if yes, 0 if no (check if PDB column is not empty)
df['has_pdb'] = df['PDB'].notna().astype(int)

# 2. Sort Data to Bubble "Best" Entries to the Top
# Primary: Gene Name (to group them)
# Secondary: Reviewed (Highest priority)
# Tertiary: Has PDB (Preferred)
# Quaternary: Length (Longest is better)
df_sorted = df.sort_values(
    by=['From', 'is_reviewed', 'has_pdb', 'Length'],
    ascending=[True, False, False, False]
)

# 3. Drop Duplicates
# Keep only the first entry for each Gene ('From')
df_unique = df_sorted.drop_duplicates(subset=['From'], keep='first')

# 4. Save
# We only need the UniProt Entry and the Gene Name
final_df = df_unique[['From', 'Entry', 'Entry Name', 'Length', 'PDB']]
final_df.to_csv(output_file, sep='\t', index=False)

print("-" * 30)
print(f"Original Hits: {len(df)}")
print(f"Unique Genes:  {df['From'].nunique()}")
print(f"Final Targets: {len(final_df)}")
print("-" * 30)
print(f"Saved cleaned list to: {output_file}")