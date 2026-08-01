import os

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
INPUT_DIR = os.path.join(BASE_PATH, "fixed_pdbs")
OUTPUT_DIR = os.path.join(BASE_PATH, "final_pdbs_for_docking")

def process_pdb(filename):
    input_path = os.path.join(INPUT_DIR, filename)
    output_path = os.path.join(OUTPUT_DIR, filename)
    
    if not os.path.exists(input_path):
        return

    fixed_lines = []
    atom_serial = 1
    
    with open(input_path, 'r') as f:
        for line in f:
            # 1. STRIP 'CONECT' RECORDS (and 'MASTER' just in case)
            if line.startswith("CONECT") or line.startswith("MASTER"):
                continue
            
            # 2. RENUMBER ATOMS (Wrap at 99,999)
            if line.startswith(("ATOM", "HETATM")):
                # Logic: MGLTools crashes if Serial > 99,999 or uses letters (A0000)
                # We wrap: 100,000 -> 1
                
                current_serial = atom_serial
                if current_serial > 99999:
                    current_serial = current_serial % 99999
                    # Handle the exact multiple edge case
                    if current_serial == 0: current_serial = 1 # or 99999, but 1 is safer for wrapping
                
                # Slicing PDB Columns strictly
                # Cols 0-6: Record Name ("ATOM  ")
                # Cols 6-11: Serial Number (Integer, Right aligned)
                # Cols 11+: Rest of the line
                
                # {:5d} ensures strict 5-character width
                new_serial_str = f"{current_serial:5d}"
                
                # Reconstruct line
                # Careful: line[:6] is chars 0-5. line[11:] is chars 11-end.
                # We inject the 5-char serial in between.
                new_line = line[:6] + new_serial_str + line[11:]
                fixed_lines.append(new_line)
                
                atom_serial += 1
            else:
                # Keep Header, Crystal, etc.
                fixed_lines.append(line)

    with open(output_path, 'w') as f:
        f.writelines(fixed_lines)
    # print(f"Processed {filename}: {atom_serial-1} atoms.")

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    files = [f for f in os.listdir(INPUT_DIR) if f.endswith(".pdb")]
    print(f"Finalizing {len(files)} PDBs (Stripping CONECT + Renumbering)...")
    
    for i, filename in enumerate(files):
        process_pdb(filename)
        if (i+1) % 50 == 0:
            print(f"Progress: {i+1}/{len(files)}")
            
    print(f"Done! Clean files are in: {OUTPUT_DIR}")
    print("Run your mgl_convert.py on THIS folder.")

if __name__ == "__main__":
    main()