#!/usr/bin/env python3

from pathlib import Path
import gzip
import sys
from Bio.PDB import MMCIFParser

def open_maybe_gz(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")

def count_residues(chain):
    residues = []
    for res in chain:
        hetflag, resseq, icode = res.id
        if hetflag == " ":
            residues.append(resseq)
    return residues

cif_path = Path(sys.argv[1])

parser = MMCIFParser(QUIET=True)
with open_maybe_gz(cif_path) as handle:
    structure = parser.get_structure("x", handle)

model = next(structure.get_models())

print("Chains in:", cif_path)
for chain in model:
    residues = count_residues(chain)
    if not residues:
        continue
    print(
        f"Chain {chain.id}: "
        f"{len(set(residues))} residues, "
        f"range {min(residues)}-{max(residues)}"
    )