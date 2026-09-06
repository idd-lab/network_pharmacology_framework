import os
import pandas as pd
import subprocess
from multiprocessing import Pool

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
RECEPTOR_DIR = os.path.join(BASE_PATH, "ready_to_dock_pdbqt")
LIGAND_FILE = os.path.join(BASE_PATH, "compound2.pdbqt")
GRID_CSV = os.path.join(BASE_PATH, "final_docking_grids.csv")
OUTPUT_DIR = os.path.join(BASE_PATH, "docking_results")
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")

# Vina Settings
VINA_EXEC = "/usr/bin/vina"
EXHAUSTIVENESS = 64
NUM_MODES = 10
NUM_CORES_PARALLEL = 10 # Number of concurrent Vina processes

def dock_target(args):
    """
    Worker function to run a single Vina docking.
    """
    uniprot_id, grid_params = args
    
    receptor_file = os.path.join(RECEPTOR_DIR, f"{uniprot_id}.pdbqt")
    output_pose = os.path.join(OUTPUT_DIR, f"{uniprot_id}_out.pdbqt")
    log_file = os.path.join(LOG_DIR, f"{uniprot_id}.log")
    
    # Check if receptor exists (some might have failed previous steps)
    if not os.path.exists(receptor_file):
        return (uniprot_id, "SKIPPED (Missing Receptor)")

    # Construct Vina Command
    # Note: --cpu 1 ensures each process only takes 1 core, 
    # allowing us to run NUM_CORES_PARALLEL processes simultaneously.
    cmd = [
        VINA_EXEC,
        "--receptor", receptor_file,
        "--ligand", LIGAND_FILE,
        "--center_x", str(grid_params['center_x']),
        "--center_y", str(grid_params['center_y']),
        "--center_z", str(grid_params['center_z']),
        "--size_x", str(grid_params['size_x']),
        "--size_y", str(grid_params['size_y']),
        "--size_z", str(grid_params['size_z']),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--out", output_pose,
        "--cpu", "1" 
    ]
    
    try:
        # Run Vina and redirect stdout to log file
        with open(log_file, "w") as f_log:
            subprocess.run(cmd, check=True, stdout=f_log, stderr=subprocess.STDOUT)
            
        return (uniprot_id, "SUCCESS")
    except subprocess.CalledProcessError:
        return (uniprot_id, "FAILED (Vina Crash)")
    except Exception as e:
        return (uniprot_id, f"FAILED ({str(e)})")

def main():
    # Setup Directories
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
    if not os.path.exists(LOG_DIR): os.makedirs(LOG_DIR)
    
    # Check Inputs
    if not os.path.exists(LIGAND_FILE):
        print("Error: Ligand file not found. Run prepare_ligand.py first.")
        return
        
    print(f"Loading Grid Config: {GRID_CSV}")
    try:
        df = pd.read_csv(GRID_CSV)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # Prepare Tasks
    tasks = []
    print(f"Queuing {len(df)} docking jobs...")
    
    for _, row in df.iterrows():
        uid = str(row['UniProt_ID'])
        # Extract grid params
        params = {
            'center_x': row['center_x'],
            'center_y': row['center_y'],
            'center_z': row['center_z'],
            'size_x': row['size_x'],
            'size_y': row['size_y'],
            'size_z': row['size_z']
        }
        tasks.append((uid, params))

    print(f"Starting Docking on {NUM_CORES_PARALLEL} cores.")
    print(f"Config: Exhaustiveness={EXHAUSTIVENESS}, Modes={NUM_MODES}")
    
    # Execution
    results = []
    with Pool(processes=NUM_CORES_PARALLEL) as pool:
        # Use imap_unordered for realtime feedback if needed, or simple map
        for i, res in enumerate(pool.imap_unordered(dock_target, tasks)):
            results.append(res)
            # Optional: Print progress every 10 completions
            if (i+1) % 10 == 0:
                print(f"Completed {i+1}/{len(tasks)} targets...")

    # Summary
    success = sum(1 for r in results if r[1] == "SUCCESS")
    failed = len(results) - success
    
    print("-" * 50)
    print("Docking Campaign Complete.")
    print(f"Success: {success}")
    print(f"Failed:  {failed}")
    print(f"Results: {OUTPUT_DIR}")
    print(f"Logs:    {LOG_DIR}")
    
    # Identify Failures
    if failed > 0:
        print("\nFailed IDs:")
        for r in results:
            if r[1] != "SUCCESS":
                print(f"{r[0]}: {r[1]}")

if __name__ == "__main__":
    main()