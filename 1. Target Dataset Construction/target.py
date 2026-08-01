import os

# --- Configuration ---
# Base directory (converted to WSL format)
base_dir = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Dataset"

# Input Files
path_a = os.path.join(base_dir, "GroupA_predicted_crude.txt")
path_b = os.path.join(base_dir, "Group_B_depmap_crude.txt")

# Output Files
out_a_cleaned = os.path.join(base_dir, "Group_A_Final.txt")      # Cleaned A
out_b_unique = os.path.join(base_dir, "Group_B_Final_Unique.txt") # B minus A
out_overlap = os.path.join(base_dir, "Tier_1_Overlap.txt")       # Intersection

def read_ids(filepath):
    """Reads a file and returns a set of unique IDs."""
    ids = set()
    if not os.path.exists(filepath):
        print(f"Error: File not found at {filepath}")
        return ids
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                clean_line = line.strip()
                # Skip empty lines or headers if they exist
                if clean_line and "uniprot" not in clean_line.lower():
                    ids.add(clean_line)
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return ids

def save_list(filepath, id_set):
    """Saves a set of IDs to a file."""
    try:
        with open(filepath, 'w') as f:
            for uid in sorted(list(id_set)):
                f.write(f"{uid}\n")
        print(f"Saved {len(id_set)} IDs to: {os.path.basename(filepath)}")
    except Exception as e:
        print(f"Error saving to {filepath}: {e}")

def main():
    print("--- Processing Lists ---")
    
    # 1. Read Inputs
    set_a = read_ids(path_a)
    set_b = read_ids(path_b)
    
    print(f"Original A count: {len(set_a)}")
    print(f"Original B count: {len(set_b)}")
    
    # 2. Logic Operations
    # Intersection: Targets in BOTH A and B
    overlap = set_a.intersection(set_b)
    
    # B Unique: Targets in B that are NOT in A
    b_unique = set_b - set_a
    
    # A Final: Usually we keep A intact (it includes the overlap), 
    # but strictly speaking, "Group A" in your final workflow is just the predicted set.
    a_final = set_a 

    # 3. Save Outputs
    print("\n--- Saving Results ---")
    save_list(out_overlap, overlap)
    save_list(out_a_cleaned, a_final)
    save_list(out_b_unique, b_unique)
    
    # 4. Final Stats
    print("\n--- Final Statistics ---")
    print(f"Tier 1 (Overlap): {len(overlap)} targets")
    print(f"Group A (Predicted): {len(a_final)} targets")
    print(f"Group B (Novel Essentials): {len(b_unique)} targets")
    
    total_unique = len(a_final) + len(b_unique)
    print(f"Total Unique Targets to Dock: {total_unique}")

if __name__ == "__main__":
    main()