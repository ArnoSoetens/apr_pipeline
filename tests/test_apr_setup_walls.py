import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import parmed as pmd
from parmed.topologyobjects import Atom

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))

import apr_setup


class FakeRestraint:
    def __init__(self):
        self.attach = {}
        self.custom_restraint_values = {}

    def initialize(self):
        pass


def wall_test_structure(near_rim_z=-1.0, far_rim_z=4.0):
    structure = pmd.Structure()
    for group_index, rim_z in enumerate((near_rim_z, far_rim_z)):
        for atom_index in range(7):
            angle = 2.0 * 3.141592653589793 * atom_index / 7
            atom = Atom(
                name=f"O{group_index * 7 + atom_index + 1}",
                atomic_number=8,
                mass=16.0,
            )
            atom.xx = 5.0 * np.cos(angle)
            atom.xy = 5.0 * np.sin(angle)
            atom.xz = rim_z
            structure.add_atom(atom, "BCD", 1)

    for name, z, residue in (("G1", 0.0, "UNL"), ("G2", 3.0, "UNL"), ("DUM", -6.0, "DM1")):
        atom = Atom(name=name, atomic_number=6, mass=12.0)
        atom.xx = atom.xy = 0.0
        atom.xz = z
        structure.add_atom(atom, residue, 1)
    return structure


def wall_test_structure_without_host():
    structure = pmd.Structure()
    for name, z, residue in (
        ("G1", 0.0, "UNL"),
        ("G2", 3.0, "UNL"),
        ("DUM", -6.0, "DM1"),
    ):
        atom = Atom(name=name, atomic_number=6, mass=12.0)
        atom.xx = atom.xy = 0.0
        atom.xz = z
        structure.add_atom(atom, residue, 1)
    return structure


class WallRestraintTests(unittest.TestCase):
    def test_guest_masks_reverses_indices_for_negative_z_backbone(self):
        guest_atoms = [
            SimpleNamespace(idx=0, xz=4.0, residue=SimpleNamespace(name="UNL")),
            SimpleNamespace(idx=1, xz=0.0, residue=SimpleNamespace(name="UNL")),
        ]
        structure = SimpleNamespace(atoms=guest_atoms)

        with patch.object(apr_setup.pmd, "load_file", return_value=structure):
            g1, g2 = apr_setup.guest_masks(0, 1)

        self.assertEqual((g1, g2), (":UNL@2", ":UNL@1"))

    def test_dummy_anchor_validation_requires_dm1_closer_to_g1(self):
        structure = pmd.Structure()
        for name, z, residue in (
            ("C1", 0.0, "UNL"),
            ("C2", 4.0, "UNL"),
            ("DUM", -6.0, "DM1"),
        ):
            atom = Atom(name=name, atomic_number=6, mass=12.0)
            atom.xx = atom.xy = 0.0
            atom.xz = z
            structure.add_atom(atom, residue, 1)

        apr_setup.validate_dummy_anchor_geometry(
            structure, ":UNL@1", ":UNL@2"
        )
        with self.assertRaisesRegex(ValueError, "DM1 is not closer to G1"):
            apr_setup.validate_dummy_anchor_geometry(
                structure, ":UNL@2", ":UNL@1"
            )

    def test_configured_distance_wall_is_created(self):
        with (
            patch.object(apr_setup, "WALLS", {
                "distance": {
                    "enabled": True,
                    "mode": "rim",
                    "rim_atom_groups": [list(range(1, 8)), list(range(8, 15))],
                    "rim_offset": 1.0,
                    "parabolic_width": 1.0,
                    "rk2": 100.0,
                    "rk3": 100.0,
                },
            }),
            patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
        ):
            walls = apr_setup.create_wall_restraints(
                wall_test_structure(), ":DM1", ":UNL@15", ":UNL@16", 3
            )

        self.assertEqual(len(walls), 1)
        self.assertEqual((walls[0].mask1, walls[0].mask2), (":DM1", ":UNL@15"))
        self.assertEqual(
            tuple(walls[0].custom_restraint_values[key] for key in ("r1", "r2", "r3", "r4")),
            (3.0, 4.0, 11.0, 12.0),
        )
        self.assertEqual(walls[0].custom_restraint_values["rk2"], 100.0)
        self.assertEqual(walls[0].custom_restraint_values["rk3"], 100.0)
        self.assertEqual(walls[0].attach["num_windows"], 3)

    def test_rim_offset_expands_flat_region_on_both_sides(self):
        with (
            patch.object(apr_setup, "WALLS", {
                "distance": {
                    "enabled": True,
                    "mode": "rim",
                    "rim_atom_groups": [list(range(1, 8)), list(range(8, 15))],
                    "rim_offset": 3.0,
                    "parabolic_width": 1.0,
                    "rk2": 100.0,
                    "rk3": 100.0,
                },
            }),
            patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
        ):
            walls = apr_setup.create_wall_restraints(
                wall_test_structure(), ":DM1", ":UNL@15", ":UNL@16", 3
            )

        self.assertEqual(
            tuple(walls[0].custom_restraint_values[key] for key in ("r1", "r2", "r3", "r4")),
            (1.0, 2.0, 13.0, 14.0),
        )

    def test_manual_r1_through_r4_are_used_directly(self):
        with (
            patch.object(apr_setup, "WALLS", {
                "distance": {
                    "enabled": True,
                    "mode": "manual",
                    "r1": 0.5,
                    "r2": 1.0,
                    "r3": 12.0,
                    "r4": 13.0,
                    "rk2": 50.0,
                    "rk3": 75.0,
                },
            }),
            patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
        ):
            walls = apr_setup.create_wall_restraints(
                wall_test_structure_without_host(), ":DM1", ":UNL@1", ":UNL@2", 3
            )

        self.assertEqual(len(walls), 1)
        self.assertEqual(
            tuple(
                walls[0].custom_restraint_values[key]
                for key in ("r1", "r2", "r3", "r4")
            ),
            (0.5, 1.0, 12.0, 13.0),
        )
        self.assertEqual(walls[0].custom_restraint_values["rk2"], 50.0)
        self.assertEqual(walls[0].custom_restraint_values["rk3"], 75.0)

    def test_wall_planes_must_bracket_both_restraint_atoms(self):
        with (
            patch.object(apr_setup, "WALLS", {
                "distance": {
                    "enabled": True,
                    "mode": "rim",
                    "rim_atom_groups": [list(range(1, 8)), list(range(8, 15))],
                    "rim_offset": 1.0,
                    "parabolic_width": 1.0,
                    "rk2": 100.0,
                    "rk3": 100.0,
                },
            }),
            patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
        ):
            with self.assertRaisesRegex(ValueError, "do not bracket"):
                apr_setup.create_wall_restraints(
                    wall_test_structure(far_rim_z=1.5),
                    ":DM1",
                    ":UNL@15",
                    ":UNL@16",
                    3,
                )

    def test_guest_within_rim_offset_beyond_rim_is_accepted(self):
        with (
            patch.object(apr_setup, "WALLS", {
                "distance": {
                    "enabled": True,
                    "mode": "rim",
                    "rim_atom_groups": [list(range(1, 8)), list(range(8, 15))],
                    "rim_offset": 1.0,
                    "parabolic_width": 1.0,
                    "rk2": 100.0,
                    "rk3": 100.0,
                },
            }),
            patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
        ):
            walls = apr_setup.create_wall_restraints(
                wall_test_structure(far_rim_z=2.5), ":DM1", ":UNL@15", ":UNL@16", 3
            )

        self.assertAlmostEqual(walls[0].custom_restraint_values["r3"], 9.5)

    def test_disabled_wall_is_omitted(self):
        with (
            patch.object(apr_setup, "WALLS", {
                "distance": {
                    "enabled": False,
                },
            }),
        ):
            walls = apr_setup.create_wall_restraints(
                None, ":DM1", ":UNL@15", ":UNL@16", 3
            )

        self.assertEqual(walls, [])


if __name__ == "__main__":
    unittest.main()
