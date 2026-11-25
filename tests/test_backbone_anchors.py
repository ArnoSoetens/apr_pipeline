"""structure.find_backbone_chain / backbone_anchors: the G1/G2 restraint anchors."""
import sys
import unittest
from pathlib import Path

import networkx as nx
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))

from structure import backbone_anchors, find_backbone_chain


def carbon_only_chain(mol):
    """The earlier anchor search (longest carbon-carbon chain), frozen here as the reference
    for guests without a heteroatom between two carbons."""
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C']
    G = nx.Graph()
    G.add_nodes_from(carbons)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in carbons and j in carbons:
            G.add_edge(i, j)
    if len(carbons) < 2:
        return carbons
    ends = [n for n in G.nodes() if G.degree(n) <= 1]
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


# guests without a heteroatom between two carbons (from the beta-CD guest lists): PFCAs, PFSAs,
# fluorotelomers, alcohols, amines, branched, cyclic and aromatic guests
CARBON_ONLY = [
    "O=C([O-])C(F)(F)C(F)(F)C(F)(F)F",
    "O=C([O-])C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "O=S(=O)([O-])C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "O=S(=O)([O-])C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "CC(C)(F)CCC(=O)[O-]",
    "CCCCS(=O)(=O)[O-]",
    "NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "CCCCCCO",
    "CCCCCCN",
    "CC(C)CC(C)(C)C",
    "OC(=O)C1CCCCC1",
    "C1CCCCC1",
    "c1ccccc1CCC(=O)O",
    "O=C(O)c1ccc(cc1)C(=O)O",
    "CCCC(=O)[O-]",
    "C(C(F)(F)F)O",
    "CC",
]


class BackboneChainTests(unittest.TestCase):
    def test_carbon_only_guests_keep_the_earlier_chain(self):
        for smiles in CARBON_ONLY:
            mol = Chem.MolFromSmiles(smiles)
            with self.subTest(smiles=smiles):
                self.assertEqual(find_backbone_chain(mol), carbon_only_chain(mol))

    def test_explicit_hydrogens_change_nothing(self):
        for smiles in CARBON_ONLY + ["CCOCCCN", "COC(CC)CCC"]:
            mol = Chem.MolFromSmiles(smiles)
            with self.subTest(smiles=smiles):
                self.assertEqual(find_backbone_chain(Chem.AddHs(mol)), find_backbone_chain(mol))

    def test_chain_ends_on_carbon_not_on_a_carboxylate_or_sulfonate_oxygen(self):
        pfba = Chem.MolFromSmiles("O=C([O-])C(F)(F)C(F)(F)C(F)(F)F")
        self.assertEqual(set(backbone_anchors(pfba)), {1, 9})
        pfbs = Chem.MolFromSmiles("O=S(=O)([O-])C(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        self.assertEqual(set(backbone_anchors(pfbs)), {4, 13})

    def test_chain_runs_through_ether_ester_and_amine_links(self):
        for smiles, ends in (("COCCO", {0, 3}),                 # was C2-C3 (1 bond)
                             ("CCNCC", {0, 4}),                 # was C0-C1
                             ("CCOCCCN", {0, 5}),               # was C3-C5
                             ("CC(=O)OCC(O)CO", {0, 7})):       # was C4-C7
            with self.subTest(smiles=smiles):
                self.assertEqual(set(backbone_anchors(Chem.MolFromSmiles(smiles))), ends)

    def test_heteroatoms_on_the_chain_stay_inside_it(self):
        chain = find_backbone_chain(Chem.MolFromSmiles("CCOCCCN"))
        self.assertEqual(chain, [0, 1, 2, 3, 4, 5])             # C C O C C C, not the N end

    def test_equal_length_chains_prefer_fewer_heteroatoms(self):
        # C0 O1 C2(C3 C4) C5 C6 C7: C0-...-C7 and C4-...-C7 are both 5 bonds; the first pair by
        # index runs through O1, the other is all carbon
        mol = Chem.MolFromSmiles("COC(CC)CCC")
        self.assertEqual(find_backbone_chain(mol), [4, 3, 2, 5, 6, 7])

    def test_rings_give_two_carbon_ends(self):
        for smiles in ("C1CCOCC1", "C1CCCCC1", "c1ccccc1"):
            mol = Chem.MolFromSmiles(smiles)
            with self.subTest(smiles=smiles):
                ends = backbone_anchors(mol)
                self.assertEqual(len(set(ends)), 2)
                self.assertTrue(all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in ends))

    def test_a_single_carbon_has_no_anchors(self):
        for smiles in ("C", "ClC(Cl)(Cl)Cl", "FC(F)(F)S(=O)(=O)[O-]"):
            with self.subTest(smiles=smiles), self.assertRaisesRegex(ValueError, "two carbon"):
                backbone_anchors(Chem.MolFromSmiles(smiles))


if __name__ == "__main__":
    unittest.main()
