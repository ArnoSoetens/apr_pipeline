"""main.py --force: what it deletes, the confirmation, and that the guard is unchanged without it.

Temporary data folders via APR_DATA only (never data/). The real redo with AmberTools is in
tests/test_setup_end_to_end.py.
"""
import csv
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from datapoints import (COMPLEX_PDB, DONE_MARKER, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, TRAJECTORY, WINDOWS,
                        force_plan, start_stage)
from setup_runner import confirm
from test_main_guard import temporary_project
from test_run import set_up_point


def simulated_point(data, name, smi=True):
    """A set-up, simulated, analyzed point with guest-only runs and some files of the user."""
    dp = set_up_point(data, name, finished=WINDOWS)
    if smi:
        (dp / f"{name}.smi").write_text("CC\n")
    for filename in ("analysis_failed.txt", "binding_free_energy.jpg", "guest_input.json",
                     "my_guest.pdb", "notes.txt"):
        (dp / filename).write_text("{}\n")
    (dp / "results.json").write_text('{"dG_bind": -1.0, "dG_sem": 0.1}\n')
    (dp / "plots").mkdir()
    (dp / "plots" / "convergence_plot.jpg").write_text("x")
    (dp / "complex").mkdir()
    (dp / "complex" / "vac.prmtop").write_text("x")
    run = dp / "windows_implicit" / "implicit_only"
    run.mkdir(parents=True)
    (run / TRAJECTORY).write_text("x")
    (run / DONE_MARKER).write_text("")
    return dp


class ForcePlanTests(unittest.TestCase):
    def test_what_force_deletes_depends_on_the_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            dp = simulated_point(temporary_project(Path(tmp)), "iso_00001")
            simulation = {"windows", "complex", "results.json", "plots", "analysis_failed.txt",
                          "binding_free_energy.jpg"}
            self.assertEqual({p.name for p in force_plan(dp, "apr_setup")}, simulation)
            for stage in ("parametrize", "placement", "structure"):
                with self.subTest(stage=stage):
                    self.assertEqual({p.name for p in force_plan(dp, stage)}, simulation | {"windows_implicit"})
            never = {f"{dp.name}.smi", "my_guest.pdb", "guest_input.json", "notes.txt", GUEST_PDB,
                     COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD}
            self.assertFalse(never & {p.name for p in force_plan(dp, "structure")})

    def test_from_the_start_means_structure_with_a_smiles_and_placement_without(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = temporary_project(Path(tmp))
            self.assertEqual(start_stage(simulated_point(data, "iso_00001")), "structure")
            self.assertEqual(start_stage(simulated_point(data, "pdb_only", smi=False)), "placement")


class FakeStdin:
    def __init__(self, tty):
        self.tty = tty

    def isatty(self):
        return self.tty


class ConfirmTests(unittest.TestCase):
    def test_only_yes_continues(self):
        for answer, expected in (("y", True), ("yes", True), ("Y\n", True), ("", False), ("n", False),
                                 ("no", False), ("sure", False)):
            with self.subTest(answer=answer):
                self.assertIs(confirm("Continue?", stdin=FakeStdin(True), ask=lambda _: answer), expected)

    def test_without_a_terminal_only_with_yes(self):
        self.assertTrue(confirm("Continue?", assume_yes=True, stdin=FakeStdin(False)))
        with self.assertRaisesRegex(SystemExit, "add --yes"):
            confirm("Continue?", stdin=FakeStdin(False))


class MainForceTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.data = temporary_project(Path(self._tmp.name))
        self.env = {**os.environ, "APR_DATA": str(self.data)}
        self.dp = simulated_point(self.data, "iso_00001")
        simulated_point(self.data, "iso_00002")                 # stays: never selected

    def tearDown(self):
        self._tmp.cleanup()

    def main(self, *args):
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / "main.py"), *args], env=self.env,
                                capture_output=True, text=True, timeout=300, stdin=subprocess.DEVNULL)
        return result.returncode, result.stdout + result.stderr

    def test_without_a_terminal_and_without_yes_nothing_is_deleted(self):
        code, out = self.main("--force", "--from", "apr_setup", "iso_00001")
        self.assertNotEqual(code, 0)
        self.assertIn(f"--force deletes the simulation output of 1 data point(s): {len(WINDOWS)} trajectories", out)
        self.assertIn("add --yes", out)
        self.assertTrue((self.dp / "windows").is_dir())
        self.assertTrue((self.dp / "results.json").is_file())

        code, out = self.main("--yes", "iso_00001")
        self.assertEqual(code, 2)
        self.assertIn("--yes only goes with --force", out)

    def test_force_deletes_just_before_the_setup_of_that_point(self):
        code, out = self.main("--force", "--from", "apr_setup", "--yes", "iso_00001")
        self.assertIn("  force: removed windows/, complex/, results.json, plots/, analysis_failed.txt, "
                      "binding_free_energy.jpg", out)
        self.assertIn("  apr_setup", out)                       # set up again (fails here: fake guest files)
        self.assertFalse((self.dp / "results.json").exists())
        self.assertFalse((self.dp / "plots").exists())
        self.assertTrue((self.dp / "windows_implicit" / "implicit_only" / DONE_MARKER).exists())   # kept
        for name in ("guest_input.json", "my_guest.pdb", "notes.txt", "iso_00001.smi"):
            self.assertTrue((self.dp / name).exists(), name)
        self.assertTrue((self.data / "iso_00002" / "results.json").exists())
        with open(self.data / "inputs" / "results.csv") as f:          # rebuilt without iso_00001
            self.assertEqual([row["smiles_folder"] for row in csv.DictReader(f)], ["iso_00002"])

    def test_without_force_the_guard_is_unchanged(self):
        code, out = self.main("--from", "apr_setup", "iso_00001")
        self.assertEqual(code, 0, out)
        self.assertIn("SKIPPED: has simulation output", out)
        self.assertIn("main.py --force redoes it", out)
        self.assertTrue((self.dp / "results.json").exists())


if __name__ == "__main__":
    unittest.main()
