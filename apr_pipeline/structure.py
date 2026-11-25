"""
smiles2pdb.py

Python re-implementation of haddocking/smiles2pdb (https://github.com/haddocking/smiles2pdb).

The original is a C++ tool built specifically to produce a small, statically
compiled, dependency-free binary for a web portal. It uses RDKit's C++ API to:
    1. Parse a SMILES string
    2. Add explicit hydrogens
    3. Embed 3D coordinates with ETKDG (random seed 42)
    4. Optimize each conformer with the UFF force field (1000 iterations)
    5. Write all conformers to a single PDB file as separate MODELs

This script does exactly the same, using RDKit's Python API directly --
no compilation, no Docker, no system dependencies beyond RDKit itself
(`pip install rdkit`).

Usage (the structure stage, see setup_runner.py):
    python structure.py                 every data point with a <name>.smi and no guest PDB yet
    python structure.py iso_00001 ...   these data points

    smiles_to_pdb(smiles, "out.pdb", num_conformers) does one SMILES from Python.
"""

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from pathlib import Path
import argparse

from apr_config import CFG
from datapoints import refuse_simulated
from guest_input import default_charge, ensure_guest_input, update_guest_input




def smiles_to_pdb(smiles_string, output_file, num_conformers=CFG["structure"]["conformers"]):
    """
    Convert a SMILES string to a PDB file with one or more 3D conformers.

    Mirrors the original C++ tool's algorithm and console messages:
    ETKDG embedding (seed=42) followed by UFF optimization (1000 iters/conf).

    Raises ValueError on any failure (invalid SMILES, embedding failure, etc.)
    """
    print(f"Parsing SMILES: {smiles_string}")

    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        raise ValueError("Failed to parse SMILES string")

    print(f"Molecule has {mol.GetNumAtoms()} atoms")

    print("Adding hydrogens...")
    mol = Chem.AddHs(mol)
    print(f"Molecule now has {mol.GetNumAtoms()} atoms")

    print(f"Generating {num_conformers} conformer(s)...")
    params = AllChem.ETKDG()
    params.randomSeed = CFG["structure"]["random_seed"]
    params.numThreads = 0  # use all available threads

    if num_conformers == 1:
        conf_id = AllChem.EmbedMolecule(mol, params)
        if conf_id < 0:
            raise ValueError("Failed to generate 3D coordinates")
        conf_ids = [conf_id]
    else:
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, num_conformers, params))
        if not conf_ids:
            raise ValueError("Failed to generate conformers")

    print(f"Successfully generated {len(conf_ids)} conformer(s)")

    print(f"Optimizing {len(conf_ids)} conformer(s) with UFF force field...")
    for i, cid in enumerate(conf_ids):
        result = AllChem.UFFOptimizeMolecule(mol, maxIters=CFG["structure"]["uff_max_iters"], confId=cid)
        if result != 0:
            print(f"Warning: UFF optimization for conformer {i + 1} returned code {result}")
    print("Geometry optimization completed")

    print(f"Writing {len(conf_ids)} conformer(s) to PDB file: {output_file}")
    with open(output_file, "w") as f:
        for i, cid in enumerate(conf_ids):
            if len(conf_ids) > 1:
                f.write(f"MODEL     {i + 1}\n")
            f.write(Chem.MolToPDBBlock(mol, confId=cid, flavor = 2))
            if len(conf_ids) > 1:
                f.write("ENDMDL\n")

    print(f"Successfully created {output_file} with {len(conf_ids)} conformer(s)")
    return mol


def backbone_graph(mol):
    """Heavy-atom graph (no H, no F) with every non-carbon dead end stripped, repeatedly, until
    all ends are carbons: heteroatoms remain only where they link carbons (ether/ester O, amine
    N, ...), so the backbone runs through them but never ends on one (carboxylate O, sulfonate
    S/O, hydroxyl O, halogens are gone). Without such a linking heteroatom this is exactly the
    carbon-carbon graph."""
    G = nx.Graph()
    G.add_nodes_from(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 9))
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in G and j in G:
            G.add_edge(i, j)
    while True:
        dead_ends = [n for n in G.nodes()
                     if G.degree(n) <= 1 and mol.GetAtomWithIdx(n).GetAtomicNum() != 6]
        if not dead_ends:
            return G
        G.remove_nodes_from(dead_ends)


def find_backbone_chain(mol):
    """Longest chain of heavy atoms (no H or F) that starts and ends on a carbon: the longest
    shortest path between two ends of backbone_graph (in a ring system without two ends:
    between any two carbons). Ties go to the chain with fewer heteroatoms, then to the
    first pair by atom index. Returns the chain's atom indices (0-based), ends first and last.

    For a guest without a heteroatom between two carbons (alkyl and perfluoroalkyl chains with
    any end groups) this gives the same chain as the earlier carbon-only search.
    """
    G = backbone_graph(mol)
    carbons = [n for n in G.nodes() if mol.GetAtomWithIdx(n).GetAtomicNum() == 6]
    if len(carbons) < 2:
        return carbons

    def heteroatoms(path):
        return sum(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 for n in path)

    best = (0, 0, None)                                  # (bonds, heteroatoms, path)

    def consider(i, j, d):
        nonlocal best
        if d < best[0] or d == 0:
            return
        path = nx.shortest_path(G, i, j)
        if d > best[0] or heteroatoms(path) < best[1]:
            best = (d, heteroatoms(path), path)

    ends = [n for n in G.nodes() if G.degree(n) <= 1]    # all carbons (dead ends were stripped)
    if len(ends) < 2:
        lengths = dict(nx.all_pairs_shortest_path_length(G))
        for i in carbons:
            for j in carbons:
                if j > i:
                    consider(i, j, lengths.get(i, {}).get(j, 0))
    else:
        for i in ends:
            for j in ends:
                if j <= i:
                    continue
                try:
                    d = nx.shortest_path_length(G, i, j)
                except nx.NetworkXNoPath:
                    continue
                consider(i, j, d)

    if best[2] is None:
        return carbons
    return best[2]

def pick_backbone_alignment_atoms(chain):

    return chain[0], chain[-1]


def backbone_anchors(mol, source="guest"):
    """G1/G2: the two ends of the longest heavy-atom chain that ends on carbons
    (find_backbone_chain; 0-based atom indices)."""
    chain = find_backbone_chain(mol)
    if len(chain) < 2:
        raise ValueError(f"guest needs at least two carbon atoms connected through heavy atoms "
                         f"other than F: {source}")
    return pick_backbone_alignment_atoms(chain)


def run(dp):
    """Structure stage: <name>.smi -> 3D guest PDB (structure.output_pdb in config.toml).

    Fills smiles and net_charge (guest_input.default_charge) in guest_input.json while still null.
    """
    refuse_simulated(dp)
    dp = Path(dp)
    smi = dp / f"{dp.name}.smi"
    if not smi.is_file():
        raise FileNotFoundError(f"no SMILES file: {smi}")
    smiles = smi.read_text().strip()
    smiles_to_pdb(smiles, str(dp / CFG["structure"]["output_pdb"]))
    values = ensure_guest_input(dp)
    if values["smiles"] is None:
        values = update_guest_input(dp, smiles=smiles)
    if values["net_charge"] is None:
        update_guest_input(dp, net_charge=default_charge(smiles))


def main(argv=None):
    from setup_runner import stage_main
    return stage_main("structure", argv)


if __name__ == "__main__":
    raise SystemExit(main())



