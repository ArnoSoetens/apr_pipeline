import sys
import unittest
from pathlib import Path

import parmed as pmd
from parmed.topologyobjects import Atom

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import simulate_guest_apr


class GuestAprRestraintTests(unittest.TestCase):
    def make_topology(self):
        structure = pmd.Structure()
        atoms = (
            ("C1", "UNL", 1, 0.0, 0.0, 0.0, 6, 12.011),
            ("C2", "UNL", 1, 0.0, 0.0, 1.5, 6, 12.011),
            ("DUM", "DM1", 2, 0.0, 0.0, -6.0, 82, 208.0),
            ("DUM", "DM2", 3, 0.0, 0.0, -9.0, 82, 208.0),
            ("DUM", "DM3", 4, 0.0, 2.2, -11.2, 82, 208.0),
        )
        for name, resname, resnum, x, y, z, atomic_number, mass in atoms:
            atom = Atom(name=name, atomic_number=atomic_number, mass=mass)
            atom.xx, atom.xy, atom.xz = x, y, z
            structure.add_atom(atom, resname, resnum)
        return structure

    def test_p000_generates_only_guest_apr_restraints(self):
        topology = self.make_topology()

        guest_restraints = simulate_guest_apr._create_p000_restraints(
            topology, ":UNL@C1", ":UNL@C2"
        )
        lines = [
            simulate_guest_apr.amber_restraint_line(restraint, "p000")
            for restraint in guest_restraints
        ]

        self.assertEqual(len(lines), 3)
        self.assertTrue(all(line for line in lines))
        self.assertIn("iat= 3,1,", lines[0])
        self.assertIn("r2=    6.00000", lines[0])
        self.assertIn("iat= 4,3,1,", lines[1])
        self.assertIn("iat= 3,1,2,", lines[2])
        self.assertEqual(
            guest_restraints[2].restraint_type,
            simulate_guest_apr.restraints.RestraintType.Angle,
        )

    def test_leap_script_loads_guest_template_without_host(self):
        leap = simulate_guest_apr._leap_script(
            "guest_with_dummy.pdb", Path("UNL.mol2"), Path("UNL.frcmod")
        )

        self.assertIn("UNL = loadmol2", leap)
        self.assertIn("guest = loadpdb", leap)
        self.assertNotIn("BCD", leap)


if __name__ == "__main__":
    unittest.main()
