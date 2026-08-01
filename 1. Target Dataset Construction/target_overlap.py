import pandas as pd
import os

# --- Configuration ---
# Define your input files with friendly names
files_map = {
    "Swiss": r"/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Dataset/SwissTargetPrediction/SwissTargetPrediction.txt",
    "ChEMBL": r"/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Dataset/ChEMBL_Target_Prediction/ChEMBL_Cleaned_UniProt.txt",
    "MTPRC": r"/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Dataset/Multitask_Target_prediction_RDKit_ChemBL/MTPRC_UniProt.txt"
}

output_file = r"/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Dataset/Group_A_Predicted_Consensus.tsv"

def read_ids_from_file(filepath):
    """
    Reads a file and extracts UniProt IDs. 
    Assumes one ID per line OR a CSV/TSV with a 'Uniprot' column.
    """
    unique_ids = set()
    try:
        # Try reading as CSV/TSV first to handle headers
        if filepath.endswith('.txt') or filepath.endswith('.tsv') or filepath.endswith('.csv'):
            try:
                # Check if it has a header like "Uniprot ID" or similar
                df = pd.read_csv(filepath, sep=None, engine='python')
                
                # Look for a column that looks like UniProt IDs
                target_col = None
                for col in df.columns:
                    if "uniprot" in col.lower():
                        target_col = col
                        break
                
                if target_col:
                    unique_ids = set(df[target_col].dropna().astype(str).str.strip())
                    return unique_ids
            except:
                pass # Fall back to raw line reading

        # Fallback: Read line by line (raw text)
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                # Clean whitespace
                clean_line = line.strip()
                # Skip header-like lines if we are reading raw
                if clean_line and "uniprot" not in clean_line.lower():
                    unique_ids.add(clean_line)
                    
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        
    return unique_ids

def main():
    print("Reading files...")
    
    # Dictionary to store results: {UniProtID: [List of Tools]}
    master_dict = {}
    
    for tool_name, path in files_map.items():
        if not os.path.exists(path):
            print(f"WARNING: File not found: {path}")
            continue
            
        print(f"  Processing {tool_name}...")
        ids = read_ids_from_file(path)
        print(f"    Found {len(ids)} unique IDs.")
        
        for uid in ids:
            if uid not in master_dict:
                master_dict[uid] = []
            master_dict[uid].append(tool_name)

    # Convert to DataFrame
    data = []
    for uid, tools in master_dict.items():
        data.append({
            'UniProt_ID': uid,
            'Count': len(tools),
            'Sources': ", ".join(tools)
        })
    
    df = pd.DataFrame(data)
    
    # Sort by Count (descending) so overlaps are at the top
    df = df.sort_values(by=['Count', 'UniProt_ID'], ascending=[False, True])
    
    # Save
    df.to_csv(output_file, sep='\t', index=False)
    
    print("-" * 40)
    print(f"Total Unique Targets Found: {len(df)}")
    print(f"Targets with Overlap (Count >= 2): {len(df[df['Count'] >= 2])}")
    print(f"Targets in ALL 3 Tools: {len(df[df['Count'] == 3])}")
    print("-" * 40)
    print(f"Saved consensus list to: {output_file}")

if __name__ == "__main__":
    main()