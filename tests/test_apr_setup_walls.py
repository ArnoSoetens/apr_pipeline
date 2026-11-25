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
    def guest_masks_for_z(self, z_values, b1_idx, b2_idx):
        guest_atoms = [
            SimpleNamespace(idx=i, xx=0.0, xy=0.0, xz=z, residue=SimpleNamespace(name="UNL"))
            for i, z in enumerate(z_values)
        ]
        structure = SimpleNamespace(atoms=guest_atoms)
        with patch.object(apr_setup.pmd, "load_file", return_value=structure):
            return apr_setup.guest_masks(b1_idx, b2_idx)

    def test_guest_masks_keep_the_given_order_when_g1_is_lower(self):
        self.assertEqual(self.guest_masks_for_z([4.0, 0.0], 1, 0), (":UNL@2", ":UNL@1"))

    def test_guest_masks_refuse_anchors_that_point_down_instead_of_swapping(self):
        with self.assertRaisesRegex(ValueError, "--from placement"):
            self.guest_masks_for_z([4.0, 0.0], 0, 1)

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

    def centroid_walls(self, radius_offset, structure=None):
        with (
            patch.object(apr_setup, "WALLS", {
                "distance": {
                    "enabled": True,
                    "mode": "rim",
                    "guest": "centroid",
                    "rim_atom_groups": [list(range(1, 8)), list(range(8, 15))],
                    "rim_offset": 2.0,                 # not used with guest = "centroid"
                    "radius_offset": radius_offset,
                    "parabolic_width": 1.0,
                    "rk2": 100.0,
                    "rk3": 100.0,
                },
            }),
            patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
        ):
            return apr_setup.create_wall_restraints(
                structure or wall_test_structure(), ":DM1", ":UNL@15", ":UNL@16", 3
            )

    def test_centroid_wall_is_a_guest_group_at_rim_plus_radius_plus_offset(self):
        # rims 5 and 10 A above DM1; guest G1/G2 at 6 and 9 A -> centroid 7.5 A, radius 1.5 A
        (wall,) = self.centroid_walls(-1.0)
        self.assertEqual((wall.mask1, wall.mask2), (":DM1", ":UNL"))
        np.testing.assert_allclose(
            [wall.custom_restraint_values[key] for key in ("r1", "r2", "r3", "r4")],
            (3.5, 4.5, 10.5, 11.5),
        )
        (wall,) = self.centroid_walls(0.0)
        np.testing.assert_allclose(
            [wall.custom_restraint_values[key] for key in ("r2", "r3")], (3.5, 11.5)
        )

    def test_centroid_starting_outside_the_wall_is_refused(self):
        # far rim 6.5 A above DM1: r3 = 6.5 + 0.5 = 7.0 < centroid 7.5 A
        with self.assertRaisesRegex(ValueError, "centroid starts outside the distance wall"):
            self.centroid_walls(-1.0, wall_test_structure(far_rim_z=0.5))

    def test_guest_centroid_counts_hydrogens_but_the_radius_only_heavy_atoms(self):
        structure = wall_test_structure()
        hydrogen = Atom(name="H1", atomic_number=1, mass=1.008)
        hydrogen.xx = hydrogen.xy = 0.0
        hydrogen.xz = 6.0
        structure.add_atom(hydrogen, "UNL", 1)
        centroid, radius = apr_setup._guest_centroid(structure, "UNL")
        np.testing.assert_allclose(centroid, (0.0, 0.0, 3.0))   # (0 + 3 + 6) / 3, unweighted
        self.assertAlmostEqual(radius, 3.0)                      # G1 at z = 0, not the H

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


def cone_test_structure():
    """Rims of radius 5 at z=-1 (5 A above DM1) and z=4 (10 A above), DM2 at z=-9."""
    structure = wall_test_structure()
    atom = Atom(name="DUM", atomic_number=6, mass=12.0)
    atom.xx = atom.xy = 0.0
    atom.xz = -9.0
    structure.add_atom(atom, "DM2", 1)
    return structure


def cone_walls(offset=0.0, enabled=True, structure=None, guest="g1"):
    walls = {
        "distance": {"rim_atom_groups": [list(range(1, 8)), list(range(8, 15))]},
        "cone": {"enabled": enabled, "angle_offset": offset, "parabolic_width": 5.0, "rk2": 100.0,
                 "guest": guest},
    }
    with (
        patch.object(apr_setup, "WALLS", walls),
        patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
    ):
        return apr_setup.create_cone_wall(
            structure or cone_test_structure(), ":DM1", ":DM2", ":UNL@15", 3
        )


class ConeWallTests(unittest.TestCase):
    def bounds(self, wall):
        return tuple(wall.custom_restraint_values[key] for key in ("r1", "r2", "r3", "r4"))

    def test_cone_passes_through_the_near_rim_oxygens(self):
        (wall,) = cone_walls()
        r2 = 180.0 - 45.0          # radius 5 at 5 A above DM1
        np.testing.assert_allclose(self.bounds(wall), (r2 - 5.0, r2, 180.0, 180.0))
        self.assertEqual((wall.mask1, wall.mask2, wall.mask3), (":DM2", ":DM1", ":UNL@15"))
        self.assertEqual(wall.custom_restraint_values["rk2"], 100.0)
        self.assertEqual(wall.attach["num_windows"], 3)

    def test_angle_offset_widens_the_cone(self):
        (wall,) = cone_walls(offset=3.0)
        np.testing.assert_allclose(self.bounds(wall)[:2], (127.0, 132.0))

    def test_rim_group_order_does_not_matter(self):
        walls = {
            "distance": {"rim_atom_groups": [list(range(8, 15)), list(range(1, 8))]},
            "cone": {"enabled": True, "angle_offset": 0.0, "parabolic_width": 5.0, "rk2": 100.0},
        }
        with (
            patch.object(apr_setup, "WALLS", walls),
            patch.object(apr_setup.restraints, "DAT_restraint", FakeRestraint),
        ):
            (wall,) = apr_setup.create_cone_wall(cone_test_structure(), ":DM1", ":DM2", ":UNL@15", 3)
        self.assertAlmostEqual(wall.custom_restraint_values["r2"], 135.0)

    def test_disabled_cone_is_omitted(self):
        self.assertEqual(cone_walls(enabled=False, structure=object()), [])

    def test_g1_outside_the_cone_is_refused(self):
        structure = cone_test_structure()
        structure.atoms[14].xx = 6.0       # G1 at z=0: 6 A off axis, 6 A above DM1 -> 45 deg out
        with self.assertRaisesRegex(ValueError, "outside the cone"):
            cone_walls(offset=-1.0, structure=structure)

    def test_centroid_cone_is_a_guest_group_with_the_same_bounds(self):
        (wall,) = cone_walls(guest="centroid")
        self.assertEqual((wall.mask1, wall.mask2, wall.mask3), (":DM2", ":DM1", ":UNL"))
        np.testing.assert_allclose(self.bounds(wall), (130.0, 135.0, 180.0, 180.0))

    def test_centroid_cone_checks_the_centroid_not_g1(self):
        structure = cone_test_structure()
        structure.atoms[15].xx = 20.0      # G2 far off axis: centroid (10, 0) at 7.5 A -> 53 deg out
        cone_walls(structure=structure)                     # G1 is on the axis: accepted
        with self.assertRaisesRegex(ValueError, "guest centroid starts outside the cone"):
            cone_walls(structure=structure, guest="centroid")


class WallConfigTests(unittest.TestCase):
    def test_guest_must_be_g1_or_centroid_for_both_walls(self):
        import copy
        import apr_config
        for wall in ("distance", "cone"):
            cfg = copy.deepcopy(apr_config.CFG)
            cfg["apr"]["restraints"]["walls"][wall]["guest"] = "com"
            with self.subTest(wall=wall), self.assertRaisesRegex(ValueError, f"{wall}.guest"):
                apr_config._validate(cfg)

    def test_centroid_rim_wall_needs_a_finite_radius_offset_of_either_sign(self):
        import copy
        import apr_config
        cfg = copy.deepcopy(apr_config.CFG)
        distance = cfg["apr"]["restraints"]["walls"]["distance"]
        distance.update(mode="rim", guest="centroid", rim_offset=-5.0,   # rim_offset unused
                        rim_atom_groups=[list(range(1, 8)), list(range(8, 15))])
        distance["radius_offset"] = -1.0
        apr_config._validate(cfg)
        distance["radius_offset"] = float("nan")
        with self.assertRaisesRegex(ValueError, "radius_offset"):
            apr_config._validate(cfg)


if __name__ == "__main__":
    unittest.main()
