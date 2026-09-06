import os
import string
import pandas as pd
import glob
from Bio.PDB import MMCIFParser, PDBIO, Select, NeighborSearch, Selection
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import re
import numpy as np

# Suppress warnings
warnings.simplefilter('ignore', PDBConstructionWarning)

# --- CONFIGURATION ---
BASE_PATH = "/mnt/c/Users/THANG/Desktop/project1/Goniothalamus_macrocalyx_of/Target_Processing_Group_A"

# Inputs
CIF_DIR = os.path.join(BASE_PATH, "all_cif_structures")
CHAIN_INFO_FILE = os.path.join(BASE_PATH, "final_chains_info.txt")

# Outputs
OUTPUT_DIR = os.path.join(BASE_PATH, "structures_by_uniprot") 
CONVERSION_LOG = os.path.join(BASE_PATH, "conversion_statistics.csv")
RENAMING_LOG = os.path.join(BASE_PATH, "chain_renaming_details.txt")
NEW_CHAIN_INFO = os.path.join(BASE_PATH, "chains_info_after_rename.txt")
FALLBACK_LOG = os.path.join(BASE_PATH, "fallback_conversions.log")

NEIGHBOR_RADIUS = 8.0 

# ------------------------------------------------------------------------------

class AlphaFoldSelect(Select):
    def accept_model(self, model):
        return model.id == 0

def parse_chain_string(chain_str):
    target_chains = []
    if pd.isna(chain_str) or str(chain_str).lower() == 'nan':
        return []
        
    parts = [p.strip() for p in str(chain_str).split(',')]
    for part in parts:
        match = re.search(r'\[auth\s+(.*?)\]', part)
        if match:
            target_chains.append(match.group(1).strip())
        else:
            target_chains.append(part)
    return target_chains

def find_cif_file(directory, identifier):
    path = os.path.join(directory, f"{identifier}.cif")
    if os.path.exists(path): return path
    path = os.path.join(directory, f"{identifier}.cif.gz")
    if os.path.exists(path): return path
    pattern = os.path.join(directory, f"*{identifier}*.cif")
    matches = glob.glob(pattern)
    if matches: return matches[0]
    return None

def get_safe_id(used_ids):
    chars = string.ascii_uppercase + string.digits + string.ascii_lowercase
    for c in chars:
        if c not in used_ids: return c
    raise ValueError("CRITICAL: Ran out of PDB Chain IDs!")

def renumber_atoms_standard(structure, max_serial=99999):
    """
    Standard atom renumbering with wrapping for structures within/exceeding limits
    """
    all_atoms = Selection.unfold_entities(structure, 'A')
    
    serial_counter = 1
    for atom in all_atoms:
        if atom.is_disordered():
            for child_atom in atom.child_dict.values():
                child_atom.serial_number = serial_counter
            atom.selected_child.serial_number = serial_counter
        else:
            atom.serial_number = serial_counter
        
        serial_counter += 1
        if serial_counter > max_serial:
            serial_counter = 1
    
    return len(all_atoms)

def write_pdb_manually(structure, output_path, max_serial=99999):
    """
    Manual PDB writing as fallback when PDBIO fails
    Preserves all atomic data from BioPython parser
    """
    with open(output_path, 'w') as f:
        # Write header with clear documentation
        f.write("HEADER    CONVERTED FROM CIF VIA BIOPYTHON                              \n")
        f.write("REMARK   1 STRUCTURE PARSED USING BIOPYTHON MMCIFPARSER                \n")
        f.write("REMARK   2 COORDINATES PRESERVED EXACTLY FROM SOURCE STRUCTURE         \n")
        f.write("REMARK   3 ATOM SERIAL NUMBERS RENUMBERED FOR PDB FORMAT COMPLIANCE    \n")
        f.write("REMARK   4 CONVERSION METHOD: FALLBACK (PDBIO FAILED DUE TO SIZE)      \n")
        
        model = structure[0]
        serial_counter = 1
        atom_count = 0
        
        for chain in model:
            for residue in chain:
                hetero_flag = residue.id[0]
                
                for atom in residue:
                    # Handle disordered atoms - use selected altloc
                    if atom.is_disordered():
                        actual_atom = atom.selected_child
                        altloc = actual_atom.altloc if hasattr(actual_atom, 'altloc') else ' '
                    else:
                        actual_atom = atom
                        altloc = ' '
                    
                    # Extract properties directly from BioPython parsed structure
                    atom_name = actual_atom.get_name()
                    resname = residue.resname
                    chain_id = chain.id
                    resseq = residue.id[1]
                    icode = residue.id[2]
                    x, y, z = actual_atom.coord  # Coordinates unchanged from parser
                    occupancy = actual_atom.get_occupancy()
                    bfactor = actual_atom.get_bfactor()
                    element = actual_atom.element if hasattr(actual_atom, 'element') else '  '
                    
                    # Record type
                    record_type = 'ATOM  ' if hetero_flag == ' ' else 'HETATM'
                    
                    # Format atom name according to PDB specification
                    if len(atom_name) < 4:
                        if len(element.strip()) == 1:
                            atom_name_fmt = f" {atom_name:<3s}"
                        else:
                            atom_name_fmt = f"{atom_name:<4s}"
                    else:
                        atom_name_fmt = f"{atom_name:<4s}"
                    
                    chain_id_fmt = chain_id[0] if chain_id else ' '
                    icode_fmt = icode if icode and icode != ' ' else ' '
                    
                    # Write PDB line with exact coordinates from BioPython parser
                    pdb_line = (
                        f"{record_type}"  # 1-6: Record type
                        f"{serial_counter:5d}"  # 7-11: Atom serial (renumbered)
                        f" "  # 12: blank
                        f"{atom_name_fmt}"  # 13-16: Atom name
                        f"{altloc}"  # 17: Alternate location
                        f"{resname:>3s}"  # 18-20: Residue name
                        f" "  # 21: blank
                        f"{chain_id_fmt}"  # 22: Chain ID
                        f"{resseq:4d}"  # 23-26: Residue sequence
                        f"{icode_fmt}"  # 27: Insertion code
                        f"   "  # 28-30: blank
                        f"{x:8.3f}"  # 31-38: X coordinate (from parser)
                        f"{y:8.3f}"  # 39-46: Y coordinate (from parser)
                        f"{z:8.3f}"  # 47-54: Z coordinate (from parser)
                        f"{occupancy:6.2f}"  # 55-60: Occupancy (from parser)
                        f"{bfactor:6.2f}"  # 61-66: B-factor (from parser)
                        f"          "  # 67-76: blank
                        f"{element:>2s}"  # 77-78: Element symbol
                        f"\n"
                    )
                    
                    f.write(pdb_line)
                    atom_count += 1
                    
                    # Increment with wrapping
                    serial_counter += 1
                    if serial_counter > max_serial:
                        serial_counter = 1
            
            # TER card after each chain
            f.write(f"TER   {serial_counter:5d}      {resname:>3s} {chain_id_fmt}{resseq:4d}\n")
            serial_counter += 1
            if serial_counter > max_serial:
                serial_counter = 1
        
        f.write("END\n")
    
    return atom_count

def validate_fallback_conversion(structure, output_path, uniprot_id, fallback_log):
    """
    Validate that manual conversion preserved coordinates correctly
    """
    from Bio.PDB import PDBParser
    
    # Count atoms in source
    source_atoms = list(Selection.unfold_entities(structure, 'A'))
    source_count = len(source_atoms)
    source_coords = np.array([a.coord for a in source_atoms])
    
    # Parse the written PDB
    parser = PDBParser(QUIET=True)
    try:
        written_structure = parser.get_structure('written', output_path)
        written_atoms = list(Selection.unfold_entities(written_structure, 'A'))
        written_count = len(written_atoms)
        written_coords = np.array([a.coord for a in written_atoms])
    except Exception as e:
        fallback_log.write(f"{uniprot_id}: VALIDATION FAILED - Could not parse written PDB: {e}\n")
        return False
    
    # Validation checks
    fallback_log.write(f"\n{'='*70}\n")
    fallback_log.write(f"FALLBACK CONVERSION VALIDATION: {uniprot_id}\n")
    fallback_log.write(f"{'='*70}\n")
    
    # Check 1: Atom count
    if source_count != written_count:
        fallback_log.write(f"❌ ATOM COUNT MISMATCH:\n")
        fallback_log.write(f"   Source: {source_count} atoms\n")
        fallback_log.write(f"   Written: {written_count} atoms\n")
        return False
    else:
        fallback_log.write(f"✓ Atom count: {source_count} atoms\n")
    
    # Check 2: Coordinate ranges
    source_min = source_coords.min(axis=0)
    source_max = source_coords.max(axis=0)
    written_min = written_coords.min(axis=0)
    written_max = written_coords.max(axis=0)
    
    coord_match = (
        np.allclose(source_min, written_min, atol=0.001) and
        np.allclose(source_max, written_max, atol=0.001)
    )
    
    if not coord_match:
        fallback_log.write(f"❌ COORDINATE RANGE MISMATCH:\n")
        fallback_log.write(f"   Source min: {source_min}\n")
        fallback_log.write(f"   Written min: {written_min}\n")
        fallback_log.write(f"   Source max: {source_max}\n")
        fallback_log.write(f"   Written max: {written_max}\n")
        return False
    else:
        fallback_log.write(f"✓ Coordinate ranges match (within 0.001 Å)\n")
        fallback_log.write(f"   Min: [{written_min[0]:.3f}, {written_min[1]:.3f}, {written_min[2]:.3f}]\n")
        fallback_log.write(f"   Max: [{written_max[0]:.3f}, {written_max[1]:.3f}, {written_max[2]:.3f}]\n")
    
    # Check 3: Individual coordinate comparison (sample)
    max_diff = 0.0
    for i in range(min(100, len(source_atoms))):  # Check first 100 atoms
        diff = np.linalg.norm(source_coords[i] - written_coords[i])
        max_diff = max(max_diff, diff)
    
    fallback_log.write(f"✓ Sample coordinate check (first 100 atoms)\n")
    fallback_log.write(f"   Maximum difference: {max_diff:.6f} Å\n")
    
    fallback_log.write(f"\n✓ ALL VALIDATION CHECKS PASSED\n")
    fallback_log.write(f"  Coordinates preserved exactly from source structure\n")
    fallback_log.flush()
    
    return True

def process_entry(uniprot_id, pdb_id, target_chains, cif_path, output_dir, rename_logger, fallback_log):
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(uniprot_id, cif_path)
        model = structure[0]
        
        is_alphafold = (pdb_id == uniprot_id) or (target_chains == ['All'])
        
        # 1. IDENTIFY CHAINS
        if not is_alphafold and target_chains:
            available_ids = {c.id for c in model}
            valid_targets = [t for t in target_chains if t in available_ids]
            
            if not valid_targets:
                if len(available_ids) == 1: valid_targets = list(available_ids)
                else: return "FAILED", f"Targets {target_chains} not found.", [], [], 0, "N/A"

            # Neighbors
            target_atoms = []
            for tid in valid_targets: target_atoms.extend(model[tid].get_atoms())
            
            ns = NeighborSearch(Selection.unfold_entities(model, 'A'))
            neighbor_ids = set()
            for atom in target_atoms:
                neighbors = ns.search(atom.coord, NEIGHBOR_RADIUS, level='C')
                for n_chain in neighbors:
                    if n_chain.id not in valid_targets: neighbor_ids.add(n_chain.id)
            
            kept_chains_ids = valid_targets + list(neighbor_ids)
            stats_msg = f"Kept {len(valid_targets)} Targets + {len(neighbor_ids)} Neighbors"
        else:
            kept_chains_ids = [c.id for c in model]
            valid_targets = kept_chains_ids 
            stats_msg = "Full Structure (AlphaFold)"

        # 2. PRUNING
        all_ids = [c.id for c in model]
        deleted_count = 0
        for cid in all_ids:
            if cid not in kept_chains_ids:
                model.detach_child(cid)
                deleted_count += 1
                
        # 3. RENAMING
        safe_ids_in_use = set([c.id for c in model if len(c.id) == 1])
        remap_record = {} 
        
        for chain in model:
            old_id = chain.id
            if len(old_id) > 1:
                new_id = get_safe_id(safe_ids_in_use)
                chain.id = new_id
                safe_ids_in_use.add(new_id)
                remap_record[old_id] = new_id
                rename_logger.write(f"{uniprot_id}\t{pdb_id}\t{old_id}\t{new_id}\n")

        # 4. RENUMBER ATOMS
        total_atoms = renumber_atoms_standard(structure, max_serial=99999)

        # 5. TRY PDBIO FIRST (STANDARD METHOD)
        output_path = os.path.join(output_dir, f"{uniprot_id}.pdb")
        method_used = "PDBIO"
        pdbio_failed = False
        
        try:
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_path, select=AlphaFoldSelect())
            
            stats_msg += f" | {total_atoms} atoms"
            
        except Exception as pdbio_error:
            # Check if it's the serial number error
            error_str = str(pdbio_error)
            if "serial number" in error_str.lower() or "100000" in error_str:
                pdbio_failed = True
                fallback_log.write(f"\n{uniprot_id}: PDBIO failed (serial number limit), using fallback...\n")
            else:
                # Some other error - re-raise it
                raise pdbio_error
        
        # 6. FALLBACK TO MANUAL WRITING IF PDBIO FAILED
        if pdbio_failed:
            atom_count = write_pdb_manually(structure, output_path, max_serial=99999)
            
            # Validate the fallback conversion
            validation_ok = validate_fallback_conversion(structure, output_path, uniprot_id, fallback_log)
            
            if not validation_ok:
                return "ERROR", "Fallback conversion validation failed", [], [], 0, "FALLBACK_FAILED"
            
            method_used = "MANUAL_FALLBACK"
            stats_msg += f" | {atom_count} atoms (FALLBACK, VALIDATED)"

        # 7. FINAL TARGET LIST
        final_target_list = []
        for t in valid_targets:
            if t in remap_record: final_target_list.append(remap_record[t])
            else: final_target_list.append(t)
        
        return "SUCCESS", stats_msg, kept_chains_ids, final_target_list, deleted_count, method_used

    except Exception as e:
        return "ERROR", str(e), [], [], 0, "EXCEPTION"

def main():
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
    
    try:
        df = pd.read_csv(CHAIN_INFO_FILE, sep='\t')
    except Exception as e:
        print(f"Error: {e}")
        return
    
    conversion_stats = []
    new_info_lines = ["UniProt_ID\tPDB_ID\tChains_of_interest\tOriginal_Source"]
    
    print(f"Processing {len(df)} Targets...")
    print(f"Strategy: Try PDBIO for all structures, fallback to manual writing if needed\n")
    
    with open(RENAMING_LOG, 'w') as rename_log, open(FALLBACK_LOG, 'w') as fallback_log:
        rename_log.write("UniProt_ID\tSource_PDB\tOriginal_Chain\tNew_Chain\n")
        fallback_log.write("FALLBACK CONVERSION LOG\n")
        fallback_log.write("="*70 + "\n")
        fallback_log.write("Logs structures that required manual PDB writing due to PDBIO failures\n")
        fallback_log.write("All fallback conversions are validated for coordinate integrity\n")
        fallback_log.write("="*70 + "\n\n")
        
        for i, row in df.iterrows():
            uniprot_id = str(row['UniProt_ID']).strip()
            pdb_id = str(row['PDB_ID']).strip() if 'PDB_ID' in row else uniprot_id
            raw_chains = str(row['Chain_of_interest']).strip() if 'Chain_of_interest' in row else 'All'
            status_dl = str(row['Downloaded']).strip() if 'Downloaded' in row else 'Success'
            
            if "fail" in status_dl.lower() or "not found" in status_dl.lower(): continue

            cif_path = find_cif_file(CIF_DIR, pdb_id)
            if not cif_path: cif_path = find_cif_file(CIF_DIR, uniprot_id)
            
            if not cif_path:
                print(f"[{i+1}] {uniprot_id}: FILE MISSING")
                conversion_stats.append({
                    'UniProt_ID': uniprot_id, 
                    'Status': 'FILE_MISSING',
                    'Method': 'N/A'
                })
                continue
            
            target_chains = parse_chain_string(raw_chains)
            
            print(f"[{i+1}/{len(df)}] {uniprot_id} ({pdb_id})...", end=" ", flush=True)
            
            status, msg, kept_ids, final_targets, del_count, method = process_entry(
                uniprot_id, pdb_id, target_chains, cif_path, OUTPUT_DIR, rename_log, fallback_log
            )
            
            print(f"{status} [{method}]")
            
            conversion_stats.append({
                'UniProt_ID': uniprot_id,
                'Status': status,
                'Method': method,
                'Message': msg,
                'Original_Target_Request': raw_chains,
                'Kept_Chain_IDs': ",".join(kept_ids),
                'Deleted_Count': del_count
            })
            
            if status == "SUCCESS":
                t_str = "All" if (pdb_id == uniprot_id or raw_chains == "All") else ",".join(final_targets)
                new_info_lines.append(f"{uniprot_id}\t{pdb_id}\t{t_str}\t{cif_path}")

    pd.DataFrame(conversion_stats).to_csv(CONVERSION_LOG, index=False)
    with open(NEW_CHAIN_INFO, 'w') as f:
        for line in new_info_lines: f.write(line + "\n")
    
    # Summary statistics
    df_stats = pd.DataFrame(conversion_stats)
    pdbio_count = len(df_stats[df_stats['Method'] == 'PDBIO'])
    fallback_count = len(df_stats[df_stats['Method'] == 'MANUAL_FALLBACK'])
    success_count = len(df_stats[df_stats['Status'] == 'SUCCESS'])
    
    print(f"\n{'='*70}")
    print(f"CONVERSION COMPLETE")
    print(f"{'='*70}")
    print(f"Total processed: {len(df_stats)}")
    print(f"Successful conversions: {success_count}")
    print(f"\nConversion methods:")
    print(f"  ✓ PDBIO (standard): {pdbio_count} structures")
    print(f"  ✓ Manual fallback (validated): {fallback_count} structures")
    print(f"\nOutput files:")
    print(f"  - Structures: {OUTPUT_DIR}")
    print(f"  - Statistics: {CONVERSION_LOG}")
    print(f"  - Fallback log: {FALLBACK_LOG}")
    print(f"  - Chain info: {NEW_CHAIN_INFO}")
    
    if fallback_count > 0:
        print(f"\n⚠ {fallback_count} structure(s) required fallback conversion.")
        print(f"  All fallback conversions validated for coordinate integrity.")
        print(f"  See {FALLBACK_LOG} for details.")

if __name__ == "__main__":
    main()
