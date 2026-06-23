
from pathlib import Path
import glob
import gzip
import pandas as pd
import numpy as np
from Bio import PDB

def parse_pdb(pdb_path):
    """Parse a PDB file."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('target', pdb_path)
    return structure

def parse_cif_gz(cif_gz_path):
    """Parse a gzipped CIF file."""
    with gzip.open(cif_gz_path, 'rt') as f:
        parser = PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure('binder', f)
    return structure

def distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D points."""
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def get_hotspot_coords(target_structure, chain_id="A"):
    """Extract hotspot atom coordinates from target protein."""
    hotspots = {
        80: ["NZ"],      # LYS80
        81: ["OG1"],     # THR81
        83: ["OE1", "OE2"]  # GLU83
    }
    
    hotspot_coords = []
    try:
        chain = target_structure[0][chain_id]
        for res_num, atom_names in hotspots.items():
            residue = chain[res_num]
            for atom_name in atom_names:
                if atom_name in residue:
                    atom = residue[atom_name]
                    hotspot_coords.append({
                        'res_num': res_num,
                        'atom_name': atom_name,
                        'coord': atom.coord
                    })
    except Exception as e:
        print(f"Error extracting hotspots: {e}")
    
    return hotspot_coords

def get_binder_atoms(binder_structure):
    """Extract all heavy atoms from binder."""
    binder_atoms = []
    try:
        for chain in binder_structure[0]:
            for residue in chain:
                for atom in residue:
                    if atom.element != 'H':
                        binder_atoms.append(atom.coord)
    except:
        pass
    return binder_atoms

def analyze_contacts(binder_cif, hotspot_coords, cutoff=4.0):
    """Analyze contacts between binder and hotspots."""
    try:
        binder_structure = parse_cif_gz(binder_cif)
        binder_atoms = get_binder_atoms(binder_structure)
        
        if not binder_atoms:
            return None
        
        n_contacts = 0
        min_distance = float('inf')
        
        # For each hotspot atom
        for hotspot in hotspot_coords:
            hotspot_coord = hotspot['coord']
            
            # Find closest binder atom
            for binder_coord in binder_atoms:
                dist = distance(hotspot_coord, binder_coord)
                
                if dist <= cutoff:
                    n_contacts += 1
                
                if dist < min_distance:
                    min_distance = dist
        
        return {
            'n_contacts': n_contacts,
            'min_distance_A': round(min_distance, 3) if min_distance != float('inf') else 999.0
        }
    
    except Exception as e:
        print(f"Error analyzing {Path(binder_cif).name}: {e}")
        return None

def main():
    # Paths
    target_pdb = Path("inputs/2N0A_dimer_cropped_clean_first.pdb")
    binder_dir = Path("work/exp_04/diffusion_out")
    output_csv = Path("work/exp_04/scores/hotspot_contacts_corrected.csv")
    
    if not target_pdb.exists():
        raise FileNotFoundError(f"Target PDB not found: {target_pdb}")
    
    # Parse target structure
    print(f"Reading target structure: {target_pdb}")
    target_structure = parse_pdb(target_pdb)
    hotspot_coords = get_hotspot_coords(target_structure, chain_id="A")
    print(f"Found {len(hotspot_coords)} hotspot atoms")
    
    # Analyze binders
    results = []
    binder_files = sorted(glob.glob(str(binder_dir / "*.cif.gz")))
    print(f"\nAnalyzing {len(binder_files)} binder structures...")
    
    for i, binder_cif in enumerate(binder_files):
        if i % 50 == 0:
            print(f"  Processed {i}/{len(binder_files)}")
        
        contact_info = analyze_contacts(binder_cif, hotspot_coords, cutoff=4.0)
        
        if contact_info:
            results.append({
                'file': Path(binder_cif).name,
                'n_contacts': contact_info['n_contacts'],
                'min_distance_A': contact_info['min_distance_A']
            })
    
    # Save results
    df = pd.DataFrame(results)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    
    print(f"\n✓ Analyzed {len(results)} binders")
    print(f"✓ Saved to: {output_csv}")
    print(f"\nTop 10 designs (best contacts):")
    print(df.nlargest(10, 'n_contacts')[['file', 'n_contacts', 'min_distance_A']])

if __name__ == "__main__":
    main()