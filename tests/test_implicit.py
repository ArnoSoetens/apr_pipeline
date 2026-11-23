"""Guest-only runs: apr_pipeline/implicit.py and scripts/run_implicit.py.

Everything runs on temporary data folders via APR_DATA (never data/). The runs use real tleap
(AmberTools) and a fake pmemd; tests/test_setup_end_to_end.py runs them with the real pmemd.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import parmed as pmd
from parmed.topologyobjects import Atom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import implicit
from datapoints import DONE_MARKER, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, TRAJECTORY
from test_main_guard import temporary_project

PFBA_ANION = "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F"

# fake pmemd: creates the files it is asked for; FAIL_IN=<text in the folder path>:<input file name>
FAKE_PMEMD = r'''#!/bin/bash
echo "$(pwd)|$*" >> "$FAKE_DIR/calls"
input=""; output=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        -i) input="$2"; shift 2 ;;
        -o) output="$2"; touch "$2"; shift 2 ;;
        -r|-x|-inf|-e|-v) touch "$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [ -n "${FAIL_IN:-}" ] && [[ "$(pwd)" == *"${FAIL_IN%%:*}"* ]] && [ "$(basename "$input")" = "${FAIL_IN#*:}" ]; then
    printf 'output\n  ERROR: something went wrong\n' > "$output"
    exit 4
fi
'''


class GuestInputFilesTests(unittest.TestCase):
    def write_guest_files(self, directory, pdb_names=("C1", "O1"), mol2_names=("C1", "O1")):
        pdb, mol2 = directory / "guest.pdb", directory / "guest.mol2"
        pdb.write_text("".join(f"HETATM{index:5d} {name:>4s} UNL A{1:4d}    "
                               "0.000   0.000   0.000  1.00  0.00           C\n"
                               for index, name in enumerate(pdb_names, start=1)))
        mol2.write_text("@<TRIPOS>MOLECULE\nUNL\n" f"{len(mol2_names)} 0 1 0 0\nSMALL\nUSER_CHARGES\n"
                        "@<TRIPOS>ATOM\n"
                        + "".join(f"{index} {name} 0.0 0.0 0.0 C.3 1 UNL 0.0\n"
                                  for index, name in enumerate(mol2_names, start=1))
                        + "@<TRIPOS>BOND\n")
        return pdb, mol2

    def test_guest_atom_names_may_be_in_different_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            implicit.validate_guest_inputs(*self.write_guest_files(Path(tmp), mol2_names=("O1", "C1")))

    def test_mismatched_atom_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "same atom names"):
                implicit.validate_guest_inputs(*self.write_guest_files(Path(tmp), mol2_names=("C1", "N1")))


class GuestRunPiecesTests(unittest.TestCase):
    def test_implicit_only_inputs_have_no_restraints_and_continue_equilibration(self):
        import run_inputs
        from string import Template
        templates = Path(run_inputs.CFG["directories"]["templates"])
        text = {name: Template((templates / template).read_text()).substitute(run_inputs.implicit_only_values())
                for name, template in run_inputs.IMPLICIT_ONLY_IN.items()}
        run = run_inputs.CFG["run"]
        self.assertIn(f"igb =           {run['igb']}", text["minimize.in"])
        self.assertIn(f"dt =            {run['dt']}", text["equilibrate.in"])
        self.assertIn(f"temp0 =         {run['temperature']}", text["production.in"])
        self.assertIn("ntr =           0", text["production.in"])
        self.assertNotIn("DISANG", text["production.in"])
        self.assertIn("ntx =           5", text["production.in"])
        self.assertIn("irest =         1", text["production.in"])

    def test_p000_has_only_the_guest_restraints(self):
        structure = pmd.Structure()
        for name, resname, resnum, x, y, z, number, mass in (
                ("C1", "UNL", 1, 0.0, 0.0, 0.0, 6, 12.011), ("C2", "UNL", 1, 0.0, 0.0, 1.5, 6, 12.011),
                ("DUM", "DM1", 2, 0.0, 0.0, -6.0, 82, 208.0), ("DUM", "DM2", 3, 0.0, 0.0, -9.0, 82, 208.0),
                ("DUM", "DM3", 4, 0.0, 2.2, -11.2, 82, 208.0)):
            atom = Atom(name=name, atomic_number=number, mass=mass)
            atom.xx, atom.xy, atom.xz = x, y, z
            structure.add_atom(atom, resname, resnum)
        from paprika.restraints import RestraintType
        from paprika.restraints.amber import amber_restraint_line
        guest_restraints = implicit._p000_restraints(structure, ":UNL@C1", ":UNL@C2")
        lines = [amber_restraint_line(r, "p000") for r in guest_restraints]
        self.assertEqual(len(lines), 3)
        self.assertIn("iat= 3,1,", lines[0])
        self.assertIn("r2=    6.00000", lines[0])
        self.assertIn("iat= 4,3,1,", lines[1])
        self.assertIn("iat= 3,1,2,", lines[2])
        self.assertEqual(guest_restraints[2].restraint_type, RestraintType.Angle)

    def test_p000_leap_script_has_the_dummies_and_no_host(self):
        leap = implicit._leap_p000("guest_with_dummy.pdb", Path("UNL.mol2"), Path("UNL.frcmod"))
        self.assertIn("UNL = loadmol2", leap)
        self.assertIn('DM1 = loadmol2 "dm1.mol2"', leap)
        self.assertNotIn("BCD", leap)

    def test_clean_plan_keeps_topologies_trajectory_and_done_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for name in ("guest.prmtop", "guest.rst7", TRAJECTORY, DONE_MARKER, "equilibration.nc", "leap.in"):
                (folder / name).write_text("x")
            (folder / "sub").mkdir()
            (folder / "sub" / "x.prmtop").write_text("x")
            removed = {p.relative_to(folder).as_posix() for p in implicit.clean_plan(folder)}
        self.assertEqual(removed, {"guest.rst7", "equilibration.nc", "leap.in", "sub/x.prmtop"})


@unittest.skipUnless(shutil.which("tleap") and shutil.which("antechamber"), "needs AmberTools")
class RunImplicitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """One parametrized guest, made once with the real setup stages."""
        cls._tmp = tempfile.TemporaryDirectory()
        data = temporary_project(Path(cls._tmp.name) / "template")
        (data / "anion").mkdir()
        (data / "anion" / "anion.smi").write_text(PFBA_ANION + "\n")
        env = {**os.environ, "APR_DATA": str(data)}
        for stage in ("structure", "placement", "parametrize"):
            result = subprocess.run([sys.executable, str(ROOT / "apr_pipeline" / f"{stage}.py"), "anion"],
                                    env=env, capture_output=True, text=True, timeout=600)
            assert result.returncode == 0, result.stdout + result.stderr
        cls.guest = data / "anion"

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        self._run_tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._run_tmp.name)
        self.data = temporary_project(self.root)
        fake = self.root / "fake"
        (fake / "bin").mkdir(parents=True)
        (fake / "bin" / "pmemd").write_text(FAKE_PMEMD)
        (fake / "bin" / "pmemd").chmod(0o755)
        self.fake = fake
        self.env = {**os.environ, "APR_DATA": str(self.data), "FAKE_DIR": str(fake),
                    "PATH": f"{fake / 'bin'}:{os.environ['PATH']}"}

    def tearDown(self):
        self._run_tmp.cleanup()

    def point(self, name):
        dp = self.data / name
        dp.mkdir()
        for filename in (GUEST_PDB, GUEST_MOL2, GUEST_FRCMOD, "guest_input.json"):
            shutil.copy(self.guest / filename, dp / filename)
        return dp

    def run_implicit(self, *args, **extra):
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / "run_implicit.py"), *args],
                                env={**self.env, **extra}, capture_output=True, text=True, timeout=600)
        return result.returncode, result.stdout + result.stderr

    def test_a_run_needs_one_of_the_two_flags(self):
        code, out = self.run_implicit()
        self.assertEqual(code, 2)
        self.assertIn("choose --implicit-only, --p000 or both", out)

    def test_both_runs_then_never_again(self):
        dp = self.point("iso_00001")
        (self.data / "raw").mkdir()                                 # not parametrized: not picked
        code, out = self.run_implicit("--implicit-only", "--p000")
        self.assertEqual(code, 0, out)
        self.assertIn("implicit finished: 2 done, 0 failed, 0 skipped", out)
        for run, top in (("implicit_only", "guest"), ("implicit_p000", "guest_apr_p000")):
            folder = dp / "windows_implicit" / run
            self.assertTrue((folder / f"{top}.prmtop").is_file(), run)
            self.assertTrue((folder / DONE_MARKER).is_file(), run)
        p000 = dp / "windows_implicit" / "implicit_p000"
        self.assertTrue((p000 / "disang.rest").read_text().startswith("&rst iat="))
        inputs = self.data / "inputs"
        calls = (self.fake / "calls").read_text()
        self.assertIn(f"-i {inputs / 'implicit_only' / 'production.in'}", calls)
        self.assertIn(f"-i {inputs / 'production_nopbc2.in'}", calls)

        code, out = self.run_implicit("--implicit-only", "--p000", "iso_00001", "raw")
        self.assertEqual(code, 0, out)
        self.assertIn("iso_00001/implicit_only SKIPPED: already done", out)
        self.assertIn("raw SKIPPED: no UNL.pdb, UNL.mol2, UNL.frcmod (parametrize first)", out)
        self.assertIn("nothing to run", out)

    def test_a_failed_run_is_recorded_and_the_others_go_on(self):
        bad, good = self.point("bad_one"), self.point("good_one")
        code, out = self.run_implicit("--implicit-only", "bad_one", "good_one",
                                      FAIL_IN="bad_one/windows_implicit/implicit_only:production.in")
        self.assertEqual(code, 1, out)
        self.assertIn("bad_one/implicit_only FAILED: production failed with exit code 4", out)
        self.assertIn("implicit finished: 1 done, 1 failed, 0 skipped", out)
        record = (bad / "windows_implicit" / "implicit_only" / "run_failed.txt").read_text()
        self.assertIn("step: production", record)
        self.assertIn("ERROR: something went wrong", record)
        self.assertFalse((bad / "windows_implicit" / "implicit_only" / DONE_MARKER).exists())
        self.assertTrue((good / "windows_implicit" / "implicit_only" / DONE_MARKER).exists())

        code, out = self.run_implicit("--implicit-only", "bad_one")        # fixed: starts again
        self.assertEqual(code, 0, out)
        self.assertFalse((bad / "windows_implicit" / "implicit_only" / "run_failed.txt").exists())

    def test_clean_cpus_and_dry_run(self):
        points = [self.point(f"iso_0000{i}") for i in (1, 2)]
        code, out = self.run_implicit("--p000", "--dry-run")
        self.assertIn("2 guest-only run(s)", out)
        self.assertFalse(any((p / "windows_implicit").exists() for p in points))
        self.assertFalse((self.data / "inputs" / "run.env").exists())

        code, out = self.run_implicit("--p000", "--cpus", "2", "--clean")
        self.assertEqual(code, 0, out)
        self.assertIn("2 at a time", out)
        for dp in points:
            kept = sorted(p.name for p in (dp / "windows_implicit" / "implicit_p000").iterdir())
            self.assertEqual(kept, sorted(["guest_apr_p000.prmtop", TRAJECTORY, DONE_MARKER]))

    def test_cleanup_handles_finished_guest_only_runs(self):
        dp = self.point("iso_00001")
        self.run_implicit("--implicit-only", "--p000")
        unfinished = dp / "windows_implicit" / "implicit_p000"
        (unfinished / DONE_MARKER).unlink()                         # pretend it did not finish
        before = sorted(p.name for p in unfinished.iterdir())
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / "cleanup.py")], env=self.env,
                                capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("cleanup iso_00001: removed", result.stdout)    # picked without results.json
        only = sorted(p.name for p in (dp / "windows_implicit" / "implicit_only").iterdir())
        self.assertEqual(only, sorted(["guest.prmtop", TRAJECTORY, DONE_MARKER]))
        self.assertEqual(sorted(p.name for p in unfinished.iterdir()), before)


if __name__ == "__main__":
    unittest.main()
