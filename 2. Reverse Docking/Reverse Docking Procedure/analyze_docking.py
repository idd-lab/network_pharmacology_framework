import os
import glob
import pandas as pd
import re
import numpy as np

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
RESULT_DIRS = [
    os.path.join(BASE_PATH, "docking_results_0", "logs"),
    os.path.join(BASE_PATH, "docking_results_1", "logs"),
    os.path.join(BASE_PATH, "docking_results_2", "logs")
]
OUTPUT_CSV = os.path.join(BASE_PATH, "docking_summary_master.csv")

def parse_vina_log(log_path):
    """
    Parses a Vina log file to find the best binding affinity (Mode 1).
    Returns: float (affinity) or None if failed.
    """
    best_affinity = None
    try:
        with open(log_path, 'r') as f:
            lines = f.readlines()
            
        start_reading = False
        for line in lines:
            # Look for the table header divider
            if "-----+------------+----------+----------" in line:
                start_reading = True
                continue
            
            if start_reading:
                # The first line after the header is Mode 1 (Best)
                # Format:    1 |      -9.5  |      0.000 |      0.000
                parts = line.split()
                if len(parts) >= 2 and parts[0] == "1":
                    best_affinity = float(parts[1])
                    break # We only care about the best mode
    except Exception as e:
        print(f"Error reading {log_path}: {e}")
        return None
        
    return best_affinity

def main():
    print("Summarizing Docking Results...")
    
    # Dictionary to store results: {UniProt_ID: [Run0_Score, Run1_Score, Run2_Score]}
    data_map = {}
    
    # 1. Collect Data
    for run_idx, log_dir in enumerate(RESULT_DIRS):
        if not os.path.exists(log_dir):
            print(f"Warning: Directory not found: {log_dir}")
            continue
            
        print(f"Processing Run {run_idx}: {log_dir}")
        log_files = glob.glob(os.path.join(log_dir, "*.log"))
        
        for log_file in log_files:
            # Extract ID from filename (e.g., "O00763.log" -> "O00763")
            filename = os.path.basename(log_file)
            uniprot_id = os.path.splitext(filename)[0]
            
            affinity = parse_vina_log(log_file)
            
            if affinity is not None:
                if uniprot_id not in data_map:
                    data_map[uniprot_id] = [np.nan] * len(RESULT_DIRS)
                
                data_map[uniprot_id][run_idx] = affinity

    # 2. Build DataFrame
    rows = []
    for uid, scores in data_map.items():
        # Clean numpy nans for calculation
        valid_scores = [s for s in scores if not np.isnan(s)]
        
        if not valid_scores:
            continue
            
        mean_score = np.mean(valid_scores)
        std_dev = np.std(valid_scores) if len(valid_scores) > 1 else 0.0
        best_score = min(valid_scores) # The absolute best among runs
        
        row = {
            "UniProt_ID": uid,
            "Run_0": scores[0] if len(scores) > 0 else None,
            "Run_1": scores[1] if len(scores) > 1 else None,
            "Run_2": scores[2] if len(scores) > 2 else None,
            "Mean_Affinity": round(mean_score, 2),
            "Std_Dev": round(std_dev, 3),
            "Best_Affinity": best_score,
            "Replicates": len(valid_scores)
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    
    # 3. Sort by Mean Affinity (Lowest/Most Negative is best)
    df = df.sort_values(by="Mean_Affinity", ascending=True)
    
    # 4. Save
    df.to_csv(OUTPUT_CSV, index=False)
    
    print("-" * 50)
    print(f"Summary Complete.")
    print(f"Total Targets: {len(df)}")
    print(f"Top 5 Hits:")
    print(df[['UniProt_ID', 'Mean_Affinity', 'Std_Dev']].head(5).to_string(index=False))
    print(f"\nSaved to: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()