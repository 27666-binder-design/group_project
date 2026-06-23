from pathlib import Path
import glob
import gzip
import pandas as pd
import numpy as np
from Bio import PDB

def parse_cif_gz(cif_gz_path):
    """Parse a gzipped CIF file and return structure."""
    with gzip.open(cif_gz_path, 'rt') as f:
        parser = PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure('protein', f)
    return structure

def get_atom_coords(structure, chain_id, residue_num, atom_names):
    """Get coordinates for specific atoms in a residue."""
    coords = []
    try:
        chain = structure[0][chain_id]
        residue = chain[residue_num]
        for atom_name in atom_names:
            if atom_name in residue:
                coords.append(residue[atom_name].coord)
    except:
        pass
    return coords

def distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D points."""
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def analyze_hotspot_contacts(cif_gz_path, hotspot_atoms):
    """Analyze which binder atoms contact hotspot atoms."""
    try:
        structure = parse_cif_gz(cif_gz_path)
        
        # Get all binder atoms (chain A in the designed complex)
        binder_atoms = []
        try:
            for chain in structure[0]:
                for residue in chain:
                    for atom in residue:
                        binder_atoms.append(atom.coord)
        except:
            return None
        
        # Count hotspot contacts
        n_contacts = 0
        min_distance = float('inf')
        
        for hotspot_atom_coord in hotspot_atoms:
            for binder_atom_coord in binder_atoms:
                dist = distance(hotspot_atom_coord, binder_atom_coord)
                n_contacts += 1 if dist <= 4.0 else 0
                min_distance = min(min_distance, dist)
        
        return {
            'n_contacts': n_contacts,
            'min_distance_A': min_distance if min_distance != float('inf') else 999.0
        }
    except Exception as e:
        return None

def main():
    cif_dir = Path("work/exp_04/diffusion_out")
    output_csv = Path("work/exp_04/scores/hotspot_contacts_summary.csv")
    
    # Target hotspots: A80, A81, A83
    hotspots = {
        80: ["NZ"],  # A80 LYS
        81: ["OG1"],  # A81 THR
        83: ["OE1", "OE2"]  # A83 GLU
    }
    
    results = []
    cif_files = sorted(glob.glob(str(cif_dir / "*.cif.gz")))
    
    print(f"Analyzing {len(cif_files)} CIF files...")
    
    for i, cif_file in enumerate(cif_files):
        if i % 50 == 0:
            print(f"  Processed {i}/{len(cif_files)}")
        
        # Get all hotspot atom coordinates
        structure = None
        try:
            structure = parse_cif_gz(cif_file)
        except:
            continue
        
        hotspot_coords = []
        try:
            chain = structure[0]["A"]
            for res_num, atom_names in hotspots.items():
                for atom_name in atom_names:
                    if res_num in chain and atom_name in chain[res_num]:
                        hotspot_coords.append(chain[res_num][atom_name].coord)
        except:
            pass
        
        if not hotspot_coords:
            continue
        
        # Analyze contacts
        contact_info = analyze_hotspot_contacts(cif_file, hotspot_coords)
        
        if contact_info:
            results.append({
                'file': Path(cif_file).name,
                'binder_chain': 'B',
                'n_contacts': contact_info['n_contacts'],
                'min_distance_A': contact_info['min_distance_A']
            })
    
    df = pd.DataFrame(results)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    
    print(f"\nAnalyzed {len(results)} designs")
    print(f"Saved to: {output_csv}")
    print(f"\nTop 10 by minimum distance:")
    print(df.nsmallest(10, 'min_distance_A')[['file', 'n_contacts', 'min_distance_A']])

if __name__ == "__main__":
    main()