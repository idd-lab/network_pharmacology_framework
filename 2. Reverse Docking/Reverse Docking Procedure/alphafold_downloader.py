import os
import requests
import time
import pandas as pd
import re

# --- Configuration ---
input_file = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A/Group_A_no_pdb.txt"
download_dir = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A/AF_structures"
summary_file = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A/AlphaFold_Download_Log.tsv"

def get_alphafold_data(uniprot_id):
    """
    Queries the AlphaFold API.
    Returns the download URL and the AlphaFold Entry ID (e.g., AF-P00533-F1).
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data and isinstance(data, list) and len(data) > 0:
                # Get the first entry (primary model)
                entry = data[0]
                cif_url = entry.get('cifUrl')
                entryId = entry.get('entryId') # This is the specific AF code
                return cif_url, entryId
        return None, None
    except Exception as e:
        print(f"    Error querying API for {uniprot_id}: {e}")
        return None, None

def download_file(url, output_path):
    """Downloads content from a URL to a local file."""
    try:
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        print(f"    Download failed: {e}")
        return False

def main():
    # 1. Setup directories
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
        print(f"Created directory: {download_dir}")

    # 2. Read UniProt IDs
    unique_ids = set()
    try:
        with open(input_file, 'r') as f:
            for line in f:
                clean_id = line.strip()
                if clean_id:
                    unique_ids.add(clean_id)
        print(f"Found {len(unique_ids)} unique UniProt IDs to process.")
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}")
        return

    # 3. Processing Loop
    summary_data = []
    
    print("-" * 60)
    for i, uniprot_id in enumerate(unique_ids, 1):
        print(f"[{i}/{len(unique_ids)}] Processing {uniprot_id}...")
        
        status = "Unknown"
        af_id = "N/A"
        filename = "N/A"
        
        # Step A: Get Data from API
        download_url, af_id_found = get_alphafold_data(uniprot_id)
        
        if download_url and af_id_found:
            af_id = af_id_found
            # Construct filename from the AF ID (Standard: AF-UniProt-F1-model_v4.cif)
            # We use the ID returned by API to be safe
            filename = f"{af_id}-model_v4.cif" 
            save_path = os.path.join(download_dir, filename)
            
            # Check if exists
            if os.path.exists(save_path):
                print(f"    Skipping (File already exists): {filename}")
                status = "Skipped (Exists)"
            else:
                print(f"    Found URL: {download_url}")
                # Step B: Download
                if download_file(download_url, save_path):
                    print(f"    Success! Saved to: {filename}")
                    status = "Success"
                else:
                    status = "Download Failed"
        else:
            print(f"    Failed: No AlphaFold entry found via API.")
            status = "No Entry Found"
            
        # Record Summary
        summary_data.append({
            'UniProt_ID': uniprot_id,
            'AlphaFold_Entry_ID': af_id,
            'Filename': filename,
            'Status': status
        })
            
        # Be nice to the server
        time.sleep(0.5)

    # 4. Save Summary File
    print("-" * 60)
    if summary_data:
        df = pd.DataFrame(summary_data)
        # Sort for cleanliness
        df = df.sort_values(by=['Status', 'UniProt_ID'])
        
        df.to_csv(summary_file, sep='\t', index=False)
        print(f"Summary saved to: {summary_file}")
        
        # Print Quick Stats
        print(f"Total Processed: {len(df)}")
        print(f"Success: {len(df[df['Status'] == 'Success'])}")
        print(f"Not Found: {len(df[df['Status'] == 'No Entry Found'])}")

if __name__ == "__main__":
    main()