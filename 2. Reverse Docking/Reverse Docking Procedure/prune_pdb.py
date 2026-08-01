import os
import pandas as pd
import re

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
STRUCTURES_DIR = os.path.join(BASE_PATH, "structures_by_uniprot")
GRID_CSV = os.path.join(BASE_PATH, "final_docking_grids.csv")
OUTPUT_DIR = os.path.join(BASE_PATH, "pruned_pdbs") # Input for Step 1

def get_pocket_chains(pocket_atm_path):
    """Identifies chains that form the pocket."""
    chains = set()
    try:
        with open(pocket_atm_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                    chains.add(line[21])
    except Exception:
        pass
    return chains

def is_in_box(x, y, z, box):
    return (box['min_x'] <= x <= box['max_x'] and
            box['min_y'] <= y <= box['max_y'] and
            box['min_z'] <= z <= box['max_z'])

def prune_pdb_smart(uniprot_id, pdb_path, pocket_chains, box, output_path):
    with open(pdb_path, 'r') as f_in, open(output_path, 'w') as f_out:
        
        # State tracking for TER logic
        # We only write a TER if the active chain being processed was kept.
        active_chain_kept = False
        
        for line in f_in:
            # --- 1. HANDLE ATOMS (Protein Chains) ---
            if line.startswith("ATOM"):
                chain = line[21]
                if chain in pocket_chains:
                    f_out.write(line)
                    active_chain_kept = True
                else:
                    active_chain_kept = False
            
            # --- 2. HANDLE TER (Termination) ---
            elif line.startswith("TER"):
                # Case A: TER has a specific Chain ID (e.g., "TER ... A")
                if len(line) > 21 and line[21] != ' ':
                    ter_chain = line[21]
                    if ter_chain in pocket_chains:
                        f_out.write(line)
                # Case B: Bare TER (Implicitly closes the current chain)
                elif active_chain_kept:
                    f_out.write(line)
                
                # After a TER, the chain is closed. Reset state.
                active_chain_kept = False
            
            # --- 3. HANDLE HETATMS (Ligands/Water) ---
            elif line.startswith("HETATM"):
                # Check Box Logic
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    if is_in_box(x, y, z, box):
                        f_out.write(line)
                except ValueError:
                    pass # Skip malformed lines

            # --- 4. OTHERS (Header, END, etc.) ---
            elif line.startswith(("CONECT", "MASTER")):
                continue # Strip these to be safe
            else:
                f_out.write(line)

def main():
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
    
    print(f"Reading Grid Config: {GRID_CSV}")
    try:
        df = pd.read_csv(GRID_CSV)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    print(f"Pruning {len(df)} structures (Smart TER Logic)...")
    
    for _, row in df.iterrows():
        uniprot_id = str(row['UniProt_ID']).strip()
        pocket_id = int(row['Pocket_ID'])
        
        # Paths
        original_pdb = os.path.join(STRUCTURES_DIR, f"{uniprot_id}.pdb")
        pocket_atm = os.path.join(STRUCTURES_DIR, f"{uniprot_id}_out", "pockets", f"pocket{pocket_id}_atm.pdb")
        output_pdb = os.path.join(OUTPUT_DIR, f"{uniprot_id}.pdb")

        if not os.path.exists(original_pdb) or not os.path.exists(pocket_atm):
            print(f"[{uniprot_id}] Missing input files. Skipped.")
            continue

        # Grid Box
        cx, cy, cz = row['center_x'], row['center_y'], row['center_z']
        sx, sy, sz = row['size_x'], row['size_y'], row['size_z']
        box = {
            'min_x': cx - sx/2, 'max_x': cx + sx/2,
            'min_y': cy - sy/2, 'max_y': cy + sy/2,
            'min_z': cz - sz/2, 'max_z': cz + sz/2
        }

        # Run Pruning
        target_chains = get_pocket_chains(pocket_atm)
        prune_pdb_smart(uniprot_id, original_pdb, target_chains, box, output_pdb)
        
    print("-" * 50)
    print(f"Pruning Complete. Output: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()