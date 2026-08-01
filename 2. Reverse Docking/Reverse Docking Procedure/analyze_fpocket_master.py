import os
import pandas as pd
import numpy as np
import glob

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
FPOCKET_DIR = os.path.join(BASE_PATH, "structures_by_uniprot")
CHAIN_INFO_FILE = os.path.join(BASE_PATH, "chains_info_after_rename.txt")
OUTPUT_CSV = os.path.join(BASE_PATH, "final_docking_grids.csv")

# Docking Box Settings
BUFFER = 2.0  # Angstroms padding

def parse_target_chains(chain_str):
    if pd.isna(chain_str): return None
    chain_str = str(chain_str).strip()
    if chain_str.lower() in ['all', 'nan']: return 'All'
    parts = chain_str.split(',')
    return set([p.strip() for p in parts if p.strip()])

def parse_fpocket_summary(info_file):
    pockets = []
    current_pocket = None
    current_score = 0.0
    try:
        with open(info_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("Pocket") and ":" in line:
                    if current_pocket is not None:
                        pockets.append((current_pocket, current_score))
                    try:
                        current_pocket = int(line.split(":")[0].replace("Pocket", "").strip())
                        current_score = 0.0 
                    except:
                        current_pocket = None
                elif "Druggability Score" in line:
                    try:
                        current_score = float(line.split(":")[-1].strip())
                    except:
                        current_score = 0.0
        if current_pocket is not None:
            pockets.append((current_pocket, current_score))
    except Exception as e:
        print(f"Error reading info file: {e}")
        return []
    pockets.sort(key=lambda x: x[1], reverse=True)
    return pockets

def validate_pocket_contact(pocket_atm_file, target_chains):
    if target_chains == 'All': return True, {'All'}
    found_chains = set()
    try:
        with open(pocket_atm_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    # Chain ID is column 21
                    chain_id = line[21].strip()
                    found_chains.add(chain_id)
    except:
        return False, set()
    common = found_chains.intersection(target_chains)
    return (len(common) > 0), found_chains

def calculate_rigorous_box_from_pqr(pqr_path):
    """
    Scientific Method:
    Box Min = min(Coordinate - Radius)
    Box Max = max(Coordinate + Radius)
    """
    min_bounds = []
    max_bounds = []
    
    try:
        with open(pqr_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    parts = line.split()
                    # Standard PQR format:
                    # Field 0: Record Name "ATOM"
                    # Field 5,6,7: X, Y, Z
                    # Field 8: Charge
                    # Field 9: Radius (Last column)
                    try:
                        # Parse Coordinates
                        x = float(parts[5])
                        y = float(parts[6])
                        z = float(parts[7])
                        
                        # Parse Radius (fpocket puts it in the last column)
                        radius = float(parts[-1])
                        
                        # Calculate extent of this sphere
                        min_bounds.append([x - radius, y - radius, z - radius])
                        max_bounds.append([x + radius, y + radius, z + radius])
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error reading PQR {pqr_path}: {e}")
        return None

    if not min_bounds: return None

    # Calculate global bounds
    min_c = np.min(min_bounds, axis=0)
    max_c = np.max(max_bounds, axis=0)
    
    center = (min_c + max_c) / 2
    size = (max_c - min_c) + (BUFFER * 2)

    return {
        'center_x': round(center[0], 3), 'center_y': round(center[1], 3), 'center_z': round(center[2], 3),
        'size_x': round(size[0], 3), 'size_y': round(size[1], 3), 'size_z': round(size[2], 3)
    }

def main():
    print(f"Reading Chain Info: {CHAIN_INFO_FILE}")
    try:
        df = pd.read_csv(CHAIN_INFO_FILE, sep='\t')
    except Exception as e:
        print(f"Error: {e}")
        return

    results = []
    print(f"Analyzing {len(df)} targets using RIGOROUS PQR method...")

    for i, row in df.iterrows():
        uniprot_id = str(row['UniProt_ID']).strip()
        pdb_id = str(row['PDB_ID']).strip()
        raw_chains = str(row['Chains_of_interest']).strip()
        target_chains = parse_target_chains(raw_chains)
        
        base_dir = os.path.join(FPOCKET_DIR, f"{uniprot_id}_out")
        info_file = os.path.join(base_dir, f"{uniprot_id}_info.txt")
        pockets_dir = os.path.join(base_dir, "pockets")
        
        if not os.path.exists(info_file):
            print(f"[{i+1}] {uniprot_id}: SKIPPING (No info file)")
            continue
            
        ranked_pockets = parse_fpocket_summary(info_file)
        selected = None
        
        for pid, score in ranked_pockets:
            pocket_atm = os.path.join(pockets_dir, f"pocket{pid}_atm.pdb")
            pocket_vert = os.path.join(pockets_dir, f"pocket{pid}_vert.pqr")
            
            if not os.path.exists(pocket_atm) or not os.path.exists(pocket_vert): continue
                
            is_valid, found_chains = validate_pocket_contact(pocket_atm, target_chains)
            
            if is_valid:
                # Use Rigorous PQR Calculation
                box = calculate_rigorous_box_from_pqr(pocket_vert)
                if box:
                    if box['size_x'] < 1.0: continue

                    selected = {
                        'UniProt_ID': uniprot_id, 
                        'PDB_ID': pdb_id, 
                        'Pocket_ID': pid,
                        'Druggability_Score': score, 
                        'Target_Chains': raw_chains,
                        'Pocket_Touching': ",".join(found_chains),
                        **box
                    }
                    print(f"[{i+1}] {uniprot_id}: Selected Pocket {pid} (Score: {score}) Box: {box['center_x']}, {box['center_y']}, {box['center_z']}")
                    break 
        
        if selected:
            results.append(selected)
        else:
            print(f"[{i+1}] {uniprot_id}: No valid pocket found.")

    if results:
        pd.DataFrame(results).to_csv(OUTPUT_CSV, index=False)
        print(f"\nDone! Saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()