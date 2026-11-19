import sys
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
import placement


def molecule_with_coordinates(coordinates):
    molecule = Chem.MolFromSmiles("CCC")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index, (x, y, z) in enumerate(coordinates):
        conformer.SetAtomPosition(index, Point3D(x, y, z))
    molecule.AddConformer(conformer)
    return molecule


class PlacementAlignmentTests(unittest.TestCase):
    def test_selects_alignment_with_least_coordinate_movement(self):
        molecule = molecule_with_coordinates(
            [(0.0, 0.0, -1.0), (0.0, 0.0, 1.0), (0.0, 2.0, 0.0)]
        )
        coordinates = molecule.GetConformer().GetPositions()
        midpoint = (coordinates[0] + coordinates[1]) / 2.0
        centered = coordinates - midpoint
        candidates = (
            placement.rotation_matrix_from_vectors([0, 0, 1], [0, 0, 1]),
            placement.rotation_matrix_from_vectors([0, 0, 1], [0, 0, -1]),
        )
        expected_movement = min(
            np.sum((centered @ rotation.T - centered) ** 2)
            for rotation in candidates
        )

        aligned, rotation, _ = placement.place_guest_in_cavity(
            molecule, 0, 1, cavity_center=[0, 0, 0], cavity_axis=[0, 0, 1]
        )
        actual_movement = np.sum((centered @ rotation.T - centered) ** 2)

        self.assertAlmostEqual(actual_movement, expected_movement)
        np.testing.assert_allclose(
            aligned.GetConformer().GetPositions(), centered @ rotation.T
        )

    def test_inverted_inputs_keep_opposite_backbone_orientations(self):
        molecule = molecule_with_coordinates(
            [(0.0, 0.0, -1.0), (0.0, 0.0, 1.0), (0.0, 2.0, 0.0)]
        )
        inverted = molecule_with_coordinates(
            [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (0.0, -2.0, 0.0)]
        )

        aligned, _, _ = placement.place_guest_in_cavity(
            molecule, 0, 1, cavity_center=[0, 0, 0], cavity_axis=[0, 0, 1]
        )
        aligned_inverted, _, _ = placement.place_guest_in_cavity(
            inverted, 0, 1, cavity_center=[0, 0, 0], cavity_axis=[0, 0, 1]
        )

        aligned_axis = (
            aligned.GetConformer().GetPositions()[1]
            - aligned.GetConformer().GetPositions()[0]
        )
        inverted_axis = (
            aligned_inverted.GetConformer().GetPositions()[1]
            - aligned_inverted.GetConformer().GetPositions()[0]
        )
        self.assertLess(np.dot(aligned_axis, inverted_axis), 0.0)


if __name__ == "__main__":
    unittest.main()
