import os
import subprocess
from multiprocessing import Pool

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
INPUT_DIR = os.path.join(BASE_PATH, "final_pdbs_for_docking")
OUTPUT_DIR = os.path.join(BASE_PATH, "ready_to_dock_pdbqt")
LOG_FILE = os.path.join(BASE_PATH, "step3_mgltools_log.txt")

# PATHS (Verify your username path matches Step 1)
MGL_PYTHON = "/home/yourusername/miniconda3/envs/mgltools_env/bin/python"
MGL_SCRIPT = "/home/yourusername/miniconda3/envs/mgltools_env/bin/prepare_receptor4.py"

NUM_CORES = 10 

def convert(filename):
    uniprot_id = filename.replace('.pdb', '')
    input_path = os.path.join(INPUT_DIR, filename)
    output_path = os.path.join(OUTPUT_DIR, f"{uniprot_id}.pdbqt")
    
    cmd = [MGL_PYTHON, MGL_SCRIPT, "-r", input_path, "-o", output_path, "-U", "nphs_lps_waters_nonstdres"]
    
    try:
        # Capture stderr to see errors if they happen
        res = subprocess.run(cmd, capture_output=True, text=True)
        
        if res.returncode == 0 and os.path.exists(output_path) and os.path.getsize(output_path) > 1000:
            return (uniprot_id, "SUCCESS", "MGLTools")
        else:
            return (uniprot_id, "FAILED", res.stderr.strip())
    except Exception as e:
        return (uniprot_id, "FAILED", str(e))

def main():
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
    
    files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.pdb')]
    print(f"Converting {len(files)} files with MGLTools...")
    
    results = []
    with Pool(processes=NUM_CORES) as pool:
        for res in pool.imap_unordered(convert, files):
            results.append(res)
            uid, stat, msg = res
            if stat == "FAILED":
                print(f"[{uid}] FAILED: {msg}")

    # Log
    with open(LOG_FILE, 'w') as f:
        f.write("UniProt_ID\tStatus\tMethod\n")
        for uid, stat, msg in sorted(results):
            f.write(f"{uid}\t{stat}\t{msg}\n")
            
    print(f"Success: {sum(1 for r in results if r[1] == 'SUCCESS')}/{len(files)}")

if __name__ == "__main__":
    main()