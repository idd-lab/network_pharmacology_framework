#!/usr/bin/env python3
"""
Script 2: PDB Downloader with Detailed Summary Stats
"""

import pandas as pd
import requests
import os
import shutil

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"
MAPPING_FILE = os.path.join(BASE_PATH, "full_length_pdb_mapping.txt")
DOWNLOAD_DIR = os.path.join(BASE_PATH, "pdb_structures") 
CHAIN_INFO_LOG = os.path.join(BASE_PATH, "final_chains_info.txt")

# Set to 'cif'
FILE_FORMAT = 'cif' 

class PDBDownloader:
    def __init__(self, mapping_file, download_dir, chain_log_file):
        self.mapping_file = mapping_file
        self.download_dir = download_dir
        self.chain_log_file = chain_log_file
        
        os.makedirs(self.download_dir, exist_ok=True)
        self.session = requests.Session()
        self.chain_info_buffer = []

        # Stats counters
        self.stats = {
            "success": 0,
            "exists": 0,
            "alphafold": 0,
            "failed": 0,
            "errors": 0
        }

    def load_mapping(self):
        try:
            self.df = pd.read_csv(self.mapping_file, sep='\t')
            return True
        except Exception as e:
            print(f"Error loading mapping: {e}")
            return False

    def get_chains_via_graphql(self, pdb_id, uniprot_id):
        query = """
        query structure($id: String!) {
          entry(entry_id: $id) {
            polymer_entities {
              rcsb_polymer_entity_container_identifiers {
                auth_asym_ids
                reference_sequence_identifiers {
                  database_accession
                  database_name
                }
              }
            }
          }
        }
        """
        url = "https://data.rcsb.org/graphql"
        try:
            response = self.session.post(url, json={'query': query, 'variables': {'id': pdb_id}}, timeout=10)
            if response.status_code != 200: return []
            data = response.json()
            if not data or 'data' not in data or 'entry' not in data['data']: return []

            entities = data['data']['entry']['polymer_entities']
            found_chains = []

            for entity in entities:
                identifiers = entity['rcsb_polymer_entity_container_identifiers']
                refs = identifiers.get('reference_sequence_identifiers', [])
                if refs:
                    for ref in refs:
                        if ref['database_name'] == 'UniProt' and ref['database_accession'] == uniprot_id:
                            auth_ids = identifiers.get('auth_asym_ids', [])
                            found_chains.extend(auth_ids)
            return sorted(list(set(found_chains)))
        except Exception:
            return []

    def download_structure(self, pdb_id):
        filename = f"{pdb_id}.{FILE_FORMAT}.gz"
        url = f"https://files.rcsb.org/download/{filename}"
        save_path = os.path.join(self.download_dir, filename)

        if os.path.exists(save_path): return "Skipped (Exists)", filename

        try:
            r = self.session.get(url, stream=True, timeout=60)
            if r.status_code == 200:
                with open(save_path, 'wb') as f:
                    r.raw.decode_content = True
                    shutil.copyfileobj(r.raw, f)
                return "Success", filename
            return "Failed", None
        except Exception:
            return "Error", None

    def run(self):
        print(f"Processing {len(self.df)} targets...")
        
        # Headers for the log file
        self.chain_info_buffer.append("UniProt_ID\tPDB_ID\tChain_of_interest\tDownloaded")

        for index, row in self.df.iterrows():
            uniprot_id = str(row['UniProt_ID']).strip()
            method = str(row['Method']).strip()
            pdb_id = str(row['PDB_ID']).strip()
            
            print(f"[{index+1}] {uniprot_id}...", end=" ", flush=True)

            # --- ALPHAFOLD ---
            if method == 'AlphaFold' or pdb_id == uniprot_id:
                print("AlphaFold.")
                self.stats["alphafold"] += 1
                self.chain_info_buffer.append(f"{uniprot_id}\t{pdb_id}\tAll\tAlphaFold")
                continue

            # --- RCSB ---
            chains = self.get_chains_via_graphql(pdb_id, uniprot_id)
            if chains:
                chains_str = ",".join(chains)
                print(f"RCSB {pdb_id} chains:{chains_str} -> ", end="")
                
                status, _ = self.download_structure(pdb_id)
                print(status)
                
                if "Success" in status:
                    self.stats["success"] += 1
                elif "Exists" in status:
                    self.stats["exists"] += 1
                else:
                    self.stats["failed"] += 1

                self.chain_info_buffer.append(f"{uniprot_id}\t{pdb_id}\t{chains_str}\t{status}")
            else:
                print(f"Chain Error (PDB {pdb_id} mismatch).")
                self.stats["errors"] += 1
                self.chain_info_buffer.append(f"{uniprot_id}\t{pdb_id}\tERROR\tFailed_Chain_ID")

        # WRITE LOG
        with open(self.chain_log_file, 'w') as f:
            for line in self.chain_info_buffer:
                f.write(line + "\n")
        
        # --- FINAL SUMMARY ---
        print("\n" + "="*40)
        print("DOWNLOAD SUMMARY")
        print("="*40)
        print(f"Total Targets Processed:  {len(self.df)}")
        print(f"AlphaFold (Skipped):      {self.stats['alphafold']}")
        print(f"RCSB Fresh Downloads:     {self.stats['success']}")
        print(f"RCSB Already Existed:     {self.stats['exists']}")
        print(f"RCSB Failed/Errors:       {self.stats['failed'] + self.stats['errors']}")
        print("="*40)
        print(f"Log saved to: {self.chain_log_file}")

if __name__ == "__main__":
    d = PDBDownloader(MAPPING_FILE, DOWNLOAD_DIR, CHAIN_INFO_LOG)
    if d.load_mapping(): d.run()