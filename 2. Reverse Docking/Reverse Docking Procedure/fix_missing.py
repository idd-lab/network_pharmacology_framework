import os
import sys
from multiprocessing import Pool
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
INPUT_DIR = os.path.join(BASE_PATH, "cleaned_pdbs")
OUTPUT_DIR = os.path.join(BASE_PATH, "fixed_pdbs")
LOG_FILE = os.path.join(BASE_PATH, "step2_fixing_log.txt")

# Parallel Settings
NUM_CORES = 10 

def process_target(filename):
    uniprot_id = filename.replace('.pdb', '')
    input_path = os.path.join(INPUT_DIR, filename)
    output_path = os.path.join(OUTPUT_DIR, filename)
    
    try:
        # Load the text-cleaned PDB
        fixer = PDBFixer(filename=input_path)
        
        # 1. Find Missing Residues (The "Map vs Reality" check)
        # Note: If the PDB lacks SEQRES headers, this simply does nothing, which is safe.
        fixer.findMissingResidues()
        
        # 2. Find Missing Atoms (Side chains that are incomplete)
        fixer.findMissingAtoms()
        
        # 3. Add Missing Atoms (Builds the loops and side chains)
        fixer.addMissingAtoms()
        
        # 4. Add Hydrogens (pH 7.4)
        # This is crucial for docking (donors/acceptors)
        fixer.addMissingHydrogens(7.4)
        
        # Save Result
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
            
        return (uniprot_id, "SUCCESS", "Fixed atoms & Hydrogens")

    except Exception as e:
        return (uniprot_id, "FAILED", str(e))

def main():
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
    
    # Get all cleaned files
    files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.pdb')]
    
    print(f"Fixing geometry for {len(files)} structures on {NUM_CORES} cores...")
    
    results = []
    with Pool(processes=NUM_CORES) as pool:
        for res in pool.imap_unordered(process_target, files):
            results.append(res)
            uid, status, msg = res
            if status == "FAILED":
                print(f"[{uid}] FAILED: {msg}")
            # Success is silent in console to keep it clean
            
    # Write Log
    success_count = 0
    with open(LOG_FILE, 'w') as log:
        log.write("UniProt_ID\tStatus\tMessage\n")
        for uid, status, msg in sorted(results):
            log.write(f"{uid}\t{status}\t{msg}\n")
            if status == "SUCCESS":
                success_count += 1
                
    print("-" * 50)
    print(f"Geometry Fix Complete.")
    print(f"Success: {success_count}/{len(files)}")
    print(f"Log saved to: {LOG_FILE}")

if __name__ == "__main__":
    # Increase recursion limit just in case P42345 (Massive) needs it for loop building
    sys.setrecursionlimit(20000)
    main()