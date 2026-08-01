import os
import re

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
INPUT_DIR = os.path.join(BASE_PATH, "pruned_pdbs")
OUTPUT_DIR = os.path.join(BASE_PATH, "cleaned_pdbs") 

# --- WHITELIST ---
KEEP_LIST = {
    # --- Ions ---
    'ZN', 'MG', 'MN', 'FE', 'CA', 'NA', 'K', 'CL', 'NI', 'CU', 'CO', 'CD', 'HG',
    # --- Biological Cofactors ---
    'HEM', 'HEC', 'FAD', 'FMN', 'NAD', 'NAP', 'NDP', 'NADH', 'NADPH',
    'ATP', 'ADP', 'GTP', 'GDP', 'AMP', 'COA', 'ACO', 'SAH', 'SAM',
    # --- Modified Amino Acids ---
    'MSE', 'PTR', 'TPO', 'SEP', 'KCX', 'LLP', 'CSO', 'PCA', 'TYS'
}

def sanitize_line(line):
    """
    Fixes formatting errors ONLY in the Occupancy/B-factor columns.
    STRICTLY ignores coordinates (Columns 0-54).
    """
    if len(line) < 54:
        return line

    # 1. Fix the '?' ghost artifact
    if '?' in line:
        line = line.replace('?', '0.00')

    # 2. SEPARATE LINE INTO ZONES
    # PDB Coordinates (X, Y, Z) end strictly at column 54.
    prefix = line[:54]  # ATOM... X... Y... Z... (DO NOT TOUCH)
    suffix = line[54:]  # Occupancy... B-factor...
    
    # 3. Fix Fused Occupancy/B-factor
    # Logic: Look for Occupancy (e.g., 1.00) fused to B-factor (e.g., 147.30)
    # Pattern: Digit.DigitDigit (Occupancy) followed immediately by Digit or Minus
    
    # Regex breakdown:
    # (\d\.\d{2})   -> Group 1: Occupancy (e.g., 1.00)
    # (-?\d+\.\d+)  -> Group 2: B-factor (Any number with a decimal)
    bfactor_pattern = re.compile(r'(\d\.\d{2})(-?\d+\.\d+)')
    
    if bfactor_pattern.search(suffix):
        suffix = bfactor_pattern.sub(r'\1 \2', suffix)

    # Reassemble
    return prefix + suffix

def process_file(filename):
    input_path = os.path.join(INPUT_DIR, filename)
    output_path = os.path.join(OUTPUT_DIR, filename)
    
    try:
        with open(input_path, 'r') as f:
            lines = f.readlines()
            
        cleaned_lines = []
        hetatm_kept = 0
        hetatm_deleted = 0
        
        for line in lines:
            # 1. Clean formatting (Only affects tail of line)
            line = sanitize_line(line)

            # 2. Filter HETATM
            if line.startswith("HETATM"):
                # Safety check for short lines
                if len(line) > 20:
                    res_name = line[17:20].strip()
                    
                    # KEEP (Zinc, Heme, etc.)
                    if res_name in KEEP_LIST:
                        cleaned_lines.append(line)
                        hetatm_kept += 1
                    else:
                        # DELETE (Water, Inhibitors, 'X' ligands)
                        hetatm_deleted += 1
                else:
                    hetatm_deleted += 1
            
            # Keep everything else (ATOM, HEADER, END, CRYST1, etc.)
            else:
                cleaned_lines.append(line)

        # Write result
        with open(output_path, 'w') as f:
            f.writelines(cleaned_lines)
            
        return f"SUCCESS: Kept {hetatm_kept}, Deleted {hetatm_deleted}"

    except Exception as e:
        return f"ERROR: {str(e)}"

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    files = [f for f in os.listdir(INPUT_DIR) if f.endswith(".pdb")]
    print(f"Cleaning {len(files)} files into '{OUTPUT_DIR}'...")

    with open("step1_cleaning_log.txt", "w") as log:
        for i, filename in enumerate(files):
            result = process_file(filename)
            # Log every single line as requested
            print(f"[{i+1}/{len(files)}] {filename} -> {result}")
            log.write(f"{filename}\t{result}\n")
            
    print("Cleaning complete.")

if __name__ == "__main__":
    main()