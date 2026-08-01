import os

# Define paths (Converted to WSL format as requested)
input_path = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Dataset/ChEMBL_Target_Prediction/ChEMBL.txt"
output_path = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Dataset/ChEMBL_Target_Prediction/ChEMBL_Cleaned.txt"

def process_chembl_file(in_file, out_file):
    print(f"Reading from: {in_file}")
    
    unique_ids = []
    # We use a set to track duplicates while keeping list for order
    seen = set()
    
    try:
        with open(in_file, 'r') as f:
            lines = f.readlines()
            
        with open(out_file, 'w') as out:
            for line in lines:
                # Remove whitespace
                clean_line = line.strip()
                if not clean_line:
                    continue
                
                # Split by pipe if it exists, otherwise it just gives the line back
                ids = clean_line.split('|')
                
                for uni_id in ids:
                    uni_id = uni_id.strip()
                    if uni_id and uni_id not in seen:
                        out.write(f"{uni_id}\n")
                        seen.add(uni_id)
                        unique_ids.append(uni_id)
                        
        print(f"Success! Processed {len(lines)} original lines into {len(unique_ids)} unique targets.")
        print(f"Saved to: {out_file}")
        
    except FileNotFoundError:
        print(f"Error: Could not find file at {in_file}")

if __name__ == "__main__":
    process_chembl_file(input_path, output_path)