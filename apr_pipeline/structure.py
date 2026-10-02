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

Usage:
    python smiles2pdb.py <SMILES_string> <output_file.pdb> [num_conformers]

    num_conformers: Number of conformers to generate (default: 1)

Example:
    python smiles2pdb.py "Cc1c(nc2ccc(F)cc2c1C(O)=O)c3ccc(cc3)c4ccccc4" 1D3G_BRE.pdb 10
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




def smiles_to_pdb(smiles_string, output_file, num_conformers=1):
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
    params.randomSeed = 42
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
        result = AllChem.UFFOptimizeMolecule(mol, maxIters=1000, confId=cid)
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


def find_carbon_backbone_chain(mol):
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C'] #find index
    G = nx.Graph()
    G.add_nodes_from(carbons) #create graph
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in carbons and j in carbons:
            G.add_edge(i, j)

    if len(carbons) < 2: #if only 2 carbons, then no chain return the indexes
        return carbons

    ends = [n for n in G.nodes() if G.degree(n) <= 1] #find the ends
    if len(ends) < 2:
        lengths = dict(nx.all_pairs_shortest_path_length(G))
        best = (0, None, None)
        for i in carbons:
            for j in carbons:
                if j <= i:
                    continue
                d = lengths.get(i, {}).get(j, 0)
                if d > best[0]:
                    best = (d, i, j)
        if best[1] is None:
            return carbons
        return nx.shortest_path(G, best[1], best[2])

    best = (0, None, None)
    for i in ends:
        for j in ends:
            if j <= i:
                continue
            try:
                d = nx.shortest_path_length(G, i, j)
            except nx.NetworkXNoPath:
                continue
            if d > best[0]:
                best = (d, i, j)

    if best[1] is None:
        return carbons
    return nx.shortest_path(G, best[1], best[2])

def pick_backbone_alignment_atoms(chain):
    
    return chain[0], chain[-1]



