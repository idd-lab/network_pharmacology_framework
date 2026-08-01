import pandas as pd
import numpy as np
import os

# --- CONFIGURATION: PATHS ---
# Adjust these if your file locations change
PATH_A = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A/docking_summary_master_A.csv"
PATH_B = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_B/docking_summary_master_B.csv"
OUTPUT_DIR = "Analysis_Results"

# --- CONFIGURATION: THE OVERLAP LIST ---
# The 15 Targets that are both Predicted (Group A) AND Essential (DepMap)
OVERLAP_IDS = [
    'O60341', 'O60885', 'P01106', 'P06493', 'P08235', 'P09237', 'P14625', 
    'P49327', 'P49336', 'P55263', 'Q02750', 'Q07869', 'Q16288', 'Q96GD4', 'Q9UBF8'
]

def load_and_prep_data():
    """Loads and labels the datasets."""
    print("1. Loading Datasets...")
    df_a = pd.read_csv(PATH_A)
    df_b = pd.read_csv(PATH_B)
    
    df_a['Source_Group'] = 'Group A (Exploratory)'
    df_b['Source_Group'] = 'Group B (Mechanism)'
    
    return df_a, df_b

def calculate_stats(df, group_name):
    """Calculates publishable statistics for a group."""
    mean_val = df['Best_Affinity'].mean()
    std_val = df['Best_Affinity'].std()
    min_val = df['Best_Affinity'].min()
    
    print(f"\n--- Statistics for {group_name} ---")
    print(f"   N (Targets):      {len(df)}")
    print(f"   Mean Affinity:    {mean_val:.3f} kcal/mol")
    print(f"   Std Deviation:    {std_val:.3f}")
    print(f"   Strongest Binder: {min_val:.3f} kcal/mol")
    
    return mean_val, std_val

def analyze_overlap(df_a, df_b):
    """
    Tier 1 Analysis: Extracts docking scores for the Consensus/DepMap Overlap list.
    Checks both groups in case of mis-assignment, prioritizes best score.
    """
    print("\n2. Analyzing Tier 1: Overlap Candidates...")
    
    # Extract rows matching the IDs
    hits_a = df_a[df_a['UniProt_ID'].isin(OVERLAP_IDS)].copy()
    hits_b = df_b[df_b['UniProt_ID'].isin(OVERLAP_IDS)].copy()
    
    # Combine (in case some are in A and some in B)
    tier1 = pd.concat([hits_a, hits_b]).drop_duplicates(subset='UniProt_ID', keep='first')
    
    # Sort by Affinity (Lowest/Most Negative is best)
    tier1 = tier1.sort_values(by='Best_Affinity')
    return tier1

def extract_top_hits(df, n=20):
    """Extracts the top N binders."""
    return df.sort_values(by='Best_Affinity').head(n)

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    # 1. Load
    df_a, df_b = load_and_prep_data()
    
    # 2. Statistics (For "Results" section text)
    mean_a, std_a = calculate_stats(df_a, "Group A")
    mean_b, std_b = calculate_stats(df_b, "Group B")
    
    # 3. Calculate Z-Scores (Standardized Binding Strength)
    # Z = (Score - Mean) / StdDev. 
    # Note: Since lower score is better, a Negative Z-score < -2.0 is significant.
    df_a['Z_Score'] = (df_a['Best_Affinity'] - mean_a) / std_a
    df_b['Z_Score'] = (df_b['Best_Affinity'] - mean_b) / std_b
    
    # 4. TIER 1: The Overlap (Consensus + Essential)
    tier1_df = analyze_overlap(df_a, df_b)
    
    # 5. TIER 2: Mechanism Drivers (Top B)
    tier2_df = extract_top_hits(df_b, n=20)
    
    # 6. TIER 3: Off-Targets / Secondary Apps (Top A)
    tier3_df = extract_top_hits(df_a, n=20)
    
    # 7. Export for Paper
    print(f"\n3. Exporting CSVs to {OUTPUT_DIR}...")
    
    cols = ['UniProt_ID', 'Best_Affinity', 'Mean_Affinity', 'Std_Dev', 'Z_Score', 'Source_Group']
    
    t1_path = os.path.join(OUTPUT_DIR, "Tier1_Overlap_Candidates.csv")
    tier1_df[cols].to_csv(t1_path, index=False)
    print(f"   -> Tier 1 Saved: {t1_path}")
    
    t2_path = os.path.join(OUTPUT_DIR, "Tier2_Mechanism_TopHits.csv")
    tier2_df[cols].to_csv(t2_path, index=False)
    print(f"   -> Tier 2 Saved: {t2_path}")
    
    t3_path = os.path.join(OUTPUT_DIR, "Tier3_Exploratory_TopHits.csv")
    tier3_df[cols].to_csv(t3_path, index=False)
    print(f"   -> Tier 3 Saved: {t3_path}")
    
    print("\nDone. Use these CSVs for your manuscript tables.")

if __name__ == "__main__":
    main()