#!/usr/bin/env python3

from pathlib import Path
import argparse
import gzip
import math
import pandas as pd
from Bio.PDB import MMCIFParser


# -------------------------
# User settings
# -------------------------

TARGET_CHAIN = "A"

# Change this if your designed binder chain has another ID.
# Common possibilities: "B", "C", or sometimes the first non-target chain.
BINDER_CHAIN = "B"

HOTSPOTS = {
    61: ["CD", "OE1", "OE2"],          # Glu61
    62: ["CG", "CD", "OE1", "NE2"],  # Gln62
}

CONTACT_CUTOFF = 4.0  # Angstroms


# -------------------------
# Helper functions
# -------------------------

def open_maybe_gz(path):
    """Open normal text file or gzipped text file."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def is_heavy_atom(atom):
    """Return True if atom is not hydrogen."""
    element = atom.element.strip().upper()
    name = atom.get_name().strip().upper()
    return element != "H" and not name.startswith("H")


def distance(atom1, atom2):
    """Euclidean distance between two Bio.PDB atoms."""
    diff = atom1.coord - atom2.coord
    return math.sqrt((diff * diff).sum())


def get_residue_by_number(chain, resnum):
    """Find standard residue by residue number."""
    for residue in chain:
        hetflag, seqid, icode = residue.id
        if hetflag == " " and seqid == resnum:
            return residue
    return None


def choose_binder_chain(model, target_chain_id):
    """
    If BINDER_CHAIN exists, use it.
    Otherwise, choose the first protein chain that is not the target chain.
    """
    if BINDER_CHAIN in model:
        return BINDER_CHAIN, model[BINDER_CHAIN]

    for chain in model:
        if chain.id != target_chain_id:
            return chain.id, chain

    raise ValueError(f"No binder chain found. Available chains: {[c.id for c in model]}")


def parse_cif(path):
    """Parse .cif or .cif.gz with Bio.PDB MMCIFParser."""
    parser = MMCIFParser(QUIET=True)
    with open_maybe_gz(path) as handle:
        structure = parser.get_structure(Path(path).stem.replace(".cif", ""), handle)
    return structure


def check_contacts(cif_file, cutoff):
    structure = parse_cif(cif_file)

    # Use first model
    model = next(structure.get_models())

    if TARGET_CHAIN not in model:
        raise ValueError(
            f"Target chain {TARGET_CHAIN} not found. "
            f"Available chains: {[c.id for c in model]}"
        )

    target_chain = model[TARGET_CHAIN]
    binder_chain_id, binder_chain = choose_binder_chain(model, TARGET_CHAIN)

    binder_atoms = [
        atom
        for residue in binder_chain
        for atom in residue
        if is_heavy_atom(atom)
    ]

    if not binder_atoms:
        raise ValueError(f"No heavy atoms found in binder chain {binder_chain_id}")

    rows = []

    for hotspot_resnum, hotspot_atoms in HOTSPOTS.items():
        residue = get_residue_by_number(target_chain, hotspot_resnum)

        if residue is None:
            rows.append({
                "file": cif_file.name,
                "target_chain": TARGET_CHAIN,
                "binder_chain": binder_chain_id,
                "target_residue": f"{TARGET_CHAIN}{hotspot_resnum}",
                "hotspot_atom": None,
                "contact": False,
                "min_distance_A": None,
                "binder_residue": None,
                "binder_atom": None,
                "note": "target residue not found",
            })
            continue

        for hotspot_atom_name in hotspot_atoms:
            if hotspot_atom_name not in residue:
                rows.append({
                    "file": cif_file.name,
                    "target_chain": TARGET_CHAIN,
                    "binder_chain": binder_chain_id,
                    "target_residue": f"{TARGET_CHAIN}{hotspot_resnum}_{residue.resname}",
                    "hotspot_atom": hotspot_atom_name,
                    "contact": False,
                    "min_distance_A": None,
                    "binder_residue": None,
                    "binder_atom": None,
                    "note": "hotspot atom not found",
                })
                continue

            hotspot_atom = residue[hotspot_atom_name]

            min_dist = float("inf")
            closest_binder_residue = None
            closest_binder_atom = None

            for atom in binder_atoms:
                d = distance(hotspot_atom, atom)
                if d < min_dist:
                    min_dist = d
                    parent_res = atom.get_parent()
                    closest_binder_residue = (
                        f"{binder_chain_id}{parent_res.id[1]}_{parent_res.resname}"
                    )
                    closest_binder_atom = atom.get_name()

            rows.append({
                "file": cif_file.name,
                "target_chain": TARGET_CHAIN,
                "binder_chain": binder_chain_id,
                "target_residue": f"{TARGET_CHAIN}{hotspot_resnum}_{residue.resname}",
                "hotspot_atom": hotspot_atom_name,
                "contact": min_dist <= cutoff,
                "min_distance_A": round(min_dist, 3),
                "binder_residue": closest_binder_residue,
                "binder_atom": closest_binder_atom,
                "note": "",
            })

    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Check whether binder chain contacts selected target hotspot atoms in CIF/CIF.GZ files."
    )
    parser.add_argument(
        "cif_dir",
        help="Directory containing .cif or .cif.gz files."
    )
    parser.add_argument(
        "--pattern",
        default="*.cif.gz",
        help='File pattern. Default: "*.cif.gz"'
    )
    parser.add_argument(
        "--output",
        default="hotspot_contacts.csv",
        help="Detailed output CSV file."
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=CONTACT_CUTOFF,
        help="Contact cutoff in Angstroms. Default: 4.0"
    )

    args = parser.parse_args()

 

    cif_dir = Path(args.cif_dir)
    cif_files = sorted(cif_dir.glob(args.pattern))

    if not cif_files:
        raise SystemExit(f"No CIF files found in {cif_dir} with pattern {args.pattern}")

    all_rows = []

    for cif_file in cif_files:
        try:
            rows = check_contacts(cif_file, args.cutoff)
            all_rows.extend(rows)
        except Exception as e:
            all_rows.append({
                "file": cif_file.name,
                "target_chain": TARGET_CHAIN,
                "binder_chain": None,
                "target_residue": None,
                "hotspot_atom": None,
                "contact": False,
                "min_distance_A": None,
                "binder_residue": None,
                "binder_atom": None,
                "note": f"ERROR: {e}",
            })

    df = pd.DataFrame(all_rows)
    df.to_csv(args.output, index=False)

    summary = (
        df.groupby("file", dropna=False)
        .agg(
            binder_chain=("binder_chain", "first"),
            n_hotspot_atoms_checked=("hotspot_atom", "count"),
            n_contacts=("contact", "sum"),
            min_distance_A=("min_distance_A", "min"),
            notes=("note", lambda x: "; ".join(sorted(set(str(v) for v in x if str(v)))))
        )
        .reset_index()
    )

    summary_file = args.output.replace(".csv", "_summary.csv")
    summary.to_csv(summary_file, index=False)

    print(f"Wrote detailed contacts to: {args.output}")
    print(f"Wrote summary to: {summary_file}")
    print()
    print(
        summary
        .sort_values(["n_contacts", "min_distance_A"], ascending=[False, True])
        .head(30)
    )


if __name__ == "__main__":
    main()