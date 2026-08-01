#!/usr/bin/env python3
"""
Script 1: Find full-length PDB structures for UniProt IDs
Reads UniProt IDs from tier_1.txt and identifies PDB entries with full sequence coverage.
Prioritizes High Coverage first, then Best (Lowest) Resolution.
"""

import requests
import time
import os

def get_uniprot_sequence_length(uniprot_id, session):
    """Get the full sequence length from UniProt using the correct API"""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    
    try:
        response = session.get(url, timeout=30)
        if response.status_code == 200:
            text = response.text
            # Parse the text format to find sequence length
            for line in text.split('\n'):
                if line.startswith('SQ   SEQUENCE'):
                    # Format: "SQ   SEQUENCE   360 AA;  40618 MW;  ..."
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            length = int(parts[2])
                            return length
                        except ValueError:
                            pass
            print(f"  Could not parse sequence length from UniProt text for {uniprot_id}")
            return None
        else:
            print(f"  Failed to fetch UniProt data for {uniprot_id}: HTTP {response.status_code}")
            return None
    except Exception as e:
        print(f"  Error fetching UniProt data for {uniprot_id}: {e}")
        return None

def get_pdb_structures_from_uniprot(uniprot_id, session):
    """Get PDB structure information from UniProt cross-references"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        response = session.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            pdb_entries = []
            
            # Extract PDB cross-references
            if 'uniProtKBCrossReferences' in data:
                for xref in data['uniProtKBCrossReferences']:
                    if xref.get('database') == 'PDB':
                        pdb_id = xref.get('id')
                        # Extract position information
                        positions = None
                        method = None
                        resolution = None
                        chains = None
                        
                        if 'properties' in xref:
                            for prop in xref['properties']:
                                key = prop.get('key')
                                value = prop.get('value')
                                if key == 'Chains':
                                    chains = value
                                elif key == 'Method':
                                    method = value
                                elif key == 'Resolution':
                                    resolution = value
                        
                        # Parse chains to get positions
                        # Format example: "A/B/C=1-360"
                        if chains:
                            if '=' in chains:
                                chain_part, pos_part = chains.split('=', 1)
                                positions = pos_part
                            else:
                                positions = None
                        
                        pdb_entries.append({
                            'pdb_id': pdb_id,
                            'positions': positions,
                            'method': method,
                            'resolution': resolution,
                            'chains': chains
                        })
            
            return pdb_entries
        else:
            print(f"  Failed to fetch PDB data for {uniprot_id}: HTTP {response.status_code}")
            return []
    except Exception as e:
        print(f"  Error fetching PDB data for {uniprot_id}: {e}")
        return []

def parse_position_range(pos_str):
    """Parse position string like '1-360' or '240-272' to get start and end"""
    if not pos_str:
        return None, None
    
    try:
        # Handle multiple ranges separated by comma
        ranges = pos_str.split(',')
        # Take the first range for simplicity
        first_range = ranges[0].strip()
        
        if '-' in first_range:
            start, end = first_range.split('-')
            return int(start.strip()), int(end.strip())
    except Exception:
        pass
    
    return None, None

def find_full_length_pdb(uniprot_id, seq_length, pdb_entries, min_coverage=0.30):
    """
    Find best PDB entry based on:
    1. Coverage (Desc) - We want MAX coverage
    2. Resolution (Asc) - We want MIN Angstroms (Lower is better)
    """
    valid_entries = []
    
    for entry in pdb_entries:
        positions = entry.get('positions')
        
        # 1. Calculate Coverage
        coverage_pct = 0
        if positions:
            start, end = parse_position_range(positions)
            if start and end:
                coverage_length = end - start + 1
                coverage_pct = coverage_length / seq_length if seq_length > 0 else 0
        
        # 2. Parse Resolution (Handle None/Strings)
        # Default to infinity if resolution is missing/invalid so it is ranked last
        res_val = float('inf') 
        try:
            if entry.get('resolution'):
                # Clean string (e.g., "1.50 A" -> "1.50")
                res_str = str(entry['resolution']).replace('A', '').strip()
                res_val = float(res_str)
        except (ValueError, TypeError):
            pass # Keep as infinity for NMR or errors

        # Store calculated metrics in the dictionary for sorting later
        entry['coverage_pct'] = coverage_pct
        entry['parsed_resolution'] = res_val
        
        # 3. Filter by minimum coverage threshold
        if coverage_pct >= min_coverage:
            valid_entries.append(entry)
    
    if not valid_entries:
        return None

    # 4. SORTING STRATEGY (The Fix)
    # Sort by:
    #   Primary Key: coverage_pct (Descending) -> We use negative sign for this
    #   Secondary Key: parsed_resolution (Ascending) -> Lower angstrom is better
    
    valid_entries.sort(key=lambda x: (-x['coverage_pct'], x['parsed_resolution']))
    
    # Debug print to show the top choice stats (optional, good for verification)
    top_choice = valid_entries[0]
    # print(f"DEBUG: Selected {top_choice['pdb_id']} | Cov: {top_choice['coverage_pct']:.2%} | Res: {top_choice['parsed_resolution']}")
    
    return top_choice

def main():
    # Configuration
    input_file = r"/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A/Group_A_Final.txt"
    output_file = r"/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A/full_length_pdb_mapping.txt"
    
    # Read UniProt IDs
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        return
    
    with open(input_file, 'r') as f:
        uniprot_ids = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(uniprot_ids)} UniProt IDs to process\n")
    
    # Setup session
    session = requests.Session()
    session.headers.update({
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    })
    
    results = []
    
    for idx, uniprot_id in enumerate(uniprot_ids, 1):
        print(f"[{idx}/{len(uniprot_ids)}] Processing {uniprot_id}")
        
        # Get sequence length
        seq_length = get_uniprot_sequence_length(uniprot_id, session)
        if not seq_length:
            print(f"  → Could not get sequence length, marking for AlphaFold")
            results.append({
                'uniprot_id': uniprot_id,
                'seq_length': 'Unknown',
                'pdb_id': uniprot_id,
                'positions': 'N/A',
                'method': 'AlphaFold',
                'resolution': 'N/A',
                'chains': 'All',
                'coverage': 'N/A'
            })
            time.sleep(0.5)
            continue
        
        print(f"  Sequence length: {seq_length} AA")
        
        # Get PDB structures
        pdb_entries = get_pdb_structures_from_uniprot(uniprot_id, session)
        if not pdb_entries:
            print(f"  → No PDB structures found, marking for AlphaFold")
            results.append({
                'uniprot_id': uniprot_id,
                'seq_length': seq_length,
                'pdb_id': uniprot_id,
                'positions': 'N/A',
                'method': 'AlphaFold',
                'resolution': 'N/A',
                'chains': 'All',
                'coverage': 'N/A'
            })
            time.sleep(0.5)
            continue
        
        print(f"  Found {len(pdb_entries)} PDB entries")
        
        # Find full-length structure (Now with Resolution Logic!)
        best_match = find_full_length_pdb(uniprot_id, seq_length, pdb_entries)
        
        if best_match:
            print(f"  → Selected PDB: {best_match['pdb_id']} (positions: {best_match['positions']}, Res: {best_match['resolution']})")
            results.append({
                'uniprot_id': uniprot_id,
                'seq_length': seq_length,
                'pdb_id': best_match['pdb_id'],
                'positions': best_match['positions'],
                'method': best_match['method'],
                'resolution': best_match['resolution'],
                'chains': best_match['chains'],
                'coverage': f"{best_match['coverage_pct']:.2%}"
            })
        else:
            print(f"  → No full-length structure found, marking for AlphaFold")
            results.append({
                'uniprot_id': uniprot_id,
                'seq_length': seq_length,
                'pdb_id': uniprot_id,
                'positions': 'N/A',
                'method': 'AlphaFold',
                'resolution': 'N/A',
                'chains': 'All',
                'coverage': 'N/A'
            })
        
        # Be respectful to servers
        time.sleep(0.5)
    
    # Write results
    with open(output_file, 'w') as f:
        f.write("UniProt_ID\tSequence_Length\tPDB_ID\tPositions\tMethod\tResolution\tChains\tCoverage\n")
        for r in results:
            f.write(f"{r['uniprot_id']}\t{r['seq_length']}\t{r['pdb_id']}\t"
                   f"{r['positions']}\t{r['method']}\t{r['resolution']}\t{r['chains']}\t{r['coverage']}\n")
    
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"Results saved to: {output_file}")
    print(f"Total processed: {len(results)}")
    print(f"With PDB structures: {len([r for r in results if r['method'] != 'AlphaFold'])}")
    print(f"Marked for AlphaFold: {len([r for r in results if r['method'] == 'AlphaFold'])}")

if __name__ == "__main__":
    main()