import pandas as pd
import numpy as np

# --- 1. CONFIGURATION ---
work_dir = "."
expr_file = work_dir + "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
dep_file = work_dir + "D2_combined_gene_dep_scores.csv"

target_genes = [
    'RFC4', 'RAD51', 'TOP2A', 'KIF11', 'ADK', 'PNP', 'NAMPT', 
    'HARS2', 'KARS1', 'WARS1', 'CYP1A2', 'CYP1A1', 'HSD11B1', 
    'AKR1C3', 'AKR1B1', 'PGR', 'CYP19A1', 'AKT1', 'ACLY', 'EZH2', 'HDAC7'
]

# Verified mapping dictionary
cell_mapping = {
    'HL-60':      {'ach': 'ACH-000002', 'ccle': 'HL60_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'},
    'MCF-7':      {'ach': 'ACH-000019', 'ccle': 'MCF7_BREAST'},
    'MDA-MB-231': {'ach': 'ACH-000768', 'ccle': 'MDAMB231_BREAST'},
    'HeLa':       {'ach': 'ACH-001086', 'ccle': 'HELA_CERVIX'},
    'A549':       {'ach': 'ACH-000681', 'ccle': 'A549_LUNG'},
    'K562':       {'ach': 'ACH-000551', 'ccle': 'K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'},
    'Jurkat':     {'ach': 'ACH-000995', 'ccle': 'JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'},
    'Raji':       {'ach': 'ACH-000654', 'ccle': 'RAJI_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE'},
    # Include potential normal proxies just in case they exist in DepMap
    'BJ':         {'ach': 'ACH-003283', 'ccle': 'BJ_SKIN'},
    'WI-38':      {'ach': 'ACH-003284', 'ccle': 'WI38_LUNG'}
}

# --- 2. EXTRACT EXPRESSION DATA ---
print("--- EXTRACTING EXPRESSION DATA ---")
target_ach_ids = [info['ach'] for info in cell_mapping.values()]
expr_chunks = []

# Chunk size 1000 to keep memory very low
for chunk in pd.read_csv(expr_file, chunksize=1000, low_memory=False):
    filtered = chunk[chunk['ModelID'].isin(target_ach_ids)]
    if not filtered.empty:
        cols_to_keep = ['ModelID']
        for col in filtered.columns:
            for gene in target_genes:
                if col.startswith(f"{gene} (") or col == gene:
                    cols_to_keep.append(col)
                    break
        
        filtered = filtered[cols_to_keep]
        filtered.columns = [col.split(' ')[0] if col != 'ModelID' else col for col in filtered.columns]
        expr_chunks.append(filtered)

if expr_chunks:
    expr_df = pd.concat(expr_chunks, ignore_index=True)
    expr_melt = expr_df.melt(id_vars=['ModelID'], var_name='Gene name', value_name='Log2(TPM+1)')
    ach_to_name = {info['ach']: name for name, info in cell_mapping.items()}
    expr_melt['Cell line'] = expr_melt['ModelID'].map(ach_to_name)
    expr_melt = expr_melt.drop(columns=['ModelID'])
    print("SUCCESS: Expression data extracted.")
else:
    print("WARNING: No expression data found for target cells.")
    expr_melt = pd.DataFrame(columns=['Gene name', 'Log2(TPM+1)', 'Cell line'])


# --- 3. EXTRACT RNAi DEPENDENCY DATA ---
print("\n--- EXTRACTING RNAi DEPENDENCY ---")
target_ccle_names = [info['ccle'] for info in cell_mapping.values()]
dep_chunks = []

for chunk in pd.read_csv(dep_file, chunksize=1000, low_memory=False):
    first_col = chunk.columns[0]
    chunk = chunk.rename(columns={first_col: 'Gene name'})
    
    # DepMap lists genes as 'RFC4 (5984)'
    chunk['Clean_Gene'] = chunk['Gene name'].apply(lambda x: str(x).split(' ')[0])
    filtered = chunk[chunk['Clean_Gene'].isin(target_genes)]
    
    if not filtered.empty:
        valid_cols = ['Clean_Gene'] + [col for col in target_ccle_names if col in filtered.columns]
        filtered = filtered[valid_cols].rename(columns={'Clean_Gene': 'Gene name'})
        dep_chunks.append(filtered)

if dep_chunks:
    dep_df = pd.concat(dep_chunks, ignore_index=True)
    dep_melt = dep_df.melt(id_vars=['Gene name'], var_name='CCLE_Name', value_name='RNAi_Dependency')
    ccle_to_name = {info['ccle']: name for name, info in cell_mapping.items()}
    dep_melt['Cell line'] = dep_melt['CCLE_Name'].map(ccle_to_name)
    dep_melt = dep_melt.drop(columns=['CCLE_Name'])
    print("SUCCESS: Dependency data extracted.")
else:
    print("WARNING: No dependency data found for target genes.")
    dep_melt = pd.DataFrame(columns=['Gene name', 'RNAi_Dependency', 'Cell line'])


# --- 4. MERGE & CALCULATE ---
print("\n--- MERGING & CALCULATING METRICS ---")
if not expr_melt.empty or not dep_melt.empty:
    final_df = pd.merge(expr_melt, dep_melt, on=['Gene name', 'Cell line'], how='outer')
    
    # We use pivot_table with aggfunc='mean' to average any duplicate gene entries
    pivot_df = final_df.pivot_table(
        index='Gene name', 
        columns='Cell line', 
        values=['Log2(TPM+1)', 'RNAi_Dependency'],
        aggfunc='mean'
    )
    
    # Flatten the hierarchical columns (e.g., ('Log2(TPM+1)', 'HL-60') -> 'HL-60_Expr')
    pivot_df.columns = [f"{cell}_{'Expr' if 'Log2' in metric else 'RNAi'}" for metric, cell in pivot_df.columns]
    
    # Identify which normal proxy actually exists in the data
    proxy = None
    for p in ['BJ', 'WI-38']:
        if f'{p}_Expr' in pivot_df.columns and not pivot_df[f'{p}_Expr'].isna().all():
            proxy = p
            break
            
    if proxy and 'HL-60_Expr' in pivot_df.columns:
        pivot_df[f'Log2FC_HL60_vs_{proxy}'] = pivot_df['HL-60_Expr'] - pivot_df[f'{proxy}_Expr']
        print(f"Calculated Expression Log2FC using {proxy} as normal baseline.")
    else:
        print("WARNING: Normal proxies (BJ/WI-38) not found. Log2FC math skipped.")

    if 'HL-60_RNAi' in pivot_df.columns and 'MCF-7_RNAi' in pivot_df.columns:
        pivot_df['Delta_RNAi_HL60_vs_MCF7'] = pivot_df['HL-60_RNAi'] - pivot_df['MCF-7_RNAi']
        print("Calculated Delta RNAi between HL-60 and MCF-7.")

    # --- 5. EXPORT ---
    output_file = work_dir + 'Q1_Unified_Omics_Matrix.csv'
    pivot_df.to_csv(output_file)
    print(f"\n--- SUCCESS ---")
    print(f"Final publication matrix saved to:\n{output_file}")
else:
    print("ERROR: Both extractions failed. Cannot merge.")