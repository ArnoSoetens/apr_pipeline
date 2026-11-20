"""Runs the real setup scripts (RDKit, antechamber, tleap, pAPRika) on a temporary project.

Never touches data/. Skipped when AmberTools is not on PATH. Takes about a minute.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from datapoints import COMPLEX_PDB, DONE, GUEST_PDB, next_stage
from test_main_guard import temporary_project

PFBA_ANION = "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F"


@unittest.skipUnless(shutil.which("tleap") and shutil.which("antechamber"), "needs AmberTools")
class SetupEndToEndTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        cls.data = temporary_project(cls.root)
        cls.env = {**os.environ, "APR_CONFIG": str(cls.root / "config.toml")}

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def run_script(self, script, *args):
        result = subprocess.run([sys.executable, str(ROOT / script), *args],
                                env=self.env, capture_output=True, text=True, timeout=600)
        return result.returncode, result.stdout + result.stderr

    def guest_input(self, name):
        return json.loads((self.data / name / "guest_input.json").read_text())

    def test_smiles_then_prepared_then_pdb_only(self):
        # 1. from a SMILES: structure -> placement -> parametrize -> apr_setup
        smiles_dp = self.data / "anion"
        smiles_dp.mkdir()
        (smiles_dp / "anion.smi").write_text(PFBA_ANION + "\n")
        code, out = self.run_script("scripts/main.py", "anion")
        self.assertEqual(code, 0, out)
        self.assertIn("1 done, 0 failed", out)
        values = self.guest_input("anion")
        self.assertEqual(values["smiles"], PFBA_ANION)
        self.assertEqual(values["net_charge"], -1)
        self.assertEqual(len(values["anchors"]), 2)
        self.assertEqual(next_stage(smiles_dp), DONE)

        # 2. prepared before guest_input.json existed: anchors are found again and saved
        prepared = self.data / "prepared"
        prepared.mkdir()
        for name in (GUEST_PDB, COMPLEX_PDB, "UNL.mol2", "UNL.frcmod"):
            shutil.copy(smiles_dp / name, prepared)
        code, out = self.run_script("scripts/main.py", "prepared")
        self.assertEqual(code, 0, out)
        self.assertIn("1 done, 0 failed", out)
        self.assertEqual(self.guest_input("prepared")["anchors"], values["anchors"])

        # 3. a PDB alone: placement.py saves anchors; parametrize.py needs a charge in the JSON
        custom = self.data / "custom"
        custom.mkdir()
        shutil.copy(smiles_dp / GUEST_PDB, custom / "my_guest.pdb")
        code, out = self.run_script("apr_pipeline/placement.py", "custom")
        self.assertEqual(code, 0, out)
        self.assertEqual(self.guest_input("custom"),
                         {"smiles": None, "net_charge": None, "anchors": values["anchors"]})

        code, out = self.run_script("apr_pipeline/parametrize.py", "custom")
        self.assertNotEqual(code, 0)
        self.assertIn("set net_charge", out)

        values_custom = self.guest_input("custom")
        values_custom["net_charge"] = -1
        (custom / "guest_input.json").write_text(json.dumps(values_custom))
        code, out = self.run_script("apr_pipeline/parametrize.py", "custom")
        self.assertEqual(code, 0, out)
        self.assertTrue((custom / "UNL.mol2").is_file())
        self.assertEqual(next_stage(custom), "apr_setup")

        # 4. a lone UNL.pdb is a placement input; any residue name in the input works
        lone = self.data / "lone"
        lone.mkdir()
        shutil.copy(smiles_dp / GUEST_PDB, lone / GUEST_PDB)
        renamed = self.data / "renamed"
        renamed.mkdir()
        (renamed / "lig.pdb").write_text((smiles_dp / GUEST_PDB).read_text().replace(" UNL ", " LIG "))
        code, out = self.run_script("apr_pipeline/placement.py", "lone", "renamed")
        self.assertEqual(code, 0, out)
        for dp in (lone, renamed):
            self.assertEqual(self.guest_input(dp.name)["anchors"], values["anchors"])
            self.assertIn(" UNL ", (dp / GUEST_PDB).read_text())
            self.assertEqual(next_stage(dp), "parametrize")

        # 5. placed before guest_input.json existed: parametrize creates it and asks for the charge
        placed = self.data / "placed"
        placed.mkdir()
        for name in (GUEST_PDB, COMPLEX_PDB):
            shutil.copy(smiles_dp / name, placed)
        code, out = self.run_script("apr_pipeline/parametrize.py", "placed")
        self.assertNotEqual(code, 0)
        self.assertIn("set net_charge", out)
        self.assertEqual(self.guest_input("placed"),
                         {"smiles": None, "net_charge": None, "anchors": None})

    def test_failure_names_the_stage_and_the_loop_goes_on(self):
        for name, smiles in (("bad", "C(("), ("good", PFBA_ANION)):
            (self.data / name).mkdir()
            (self.data / name / f"{name}.smi").write_text(smiles + "\n")
        code, out = self.run_script("scripts/main.py", "bad", "good")
        self.assertEqual(code, 0, out)
        self.assertIn("bad FAILED in structure:", out)
        self.assertIn("1 done, 1 failed", out)
        self.assertTrue((self.data / "bad" / "setup_failed.txt").read_text()
                        .splitlines()[1].startswith("structure: "))
        self.assertFalse((self.data / "good" / "setup_failed.txt").exists())


if __name__ == "__main__":
    unittest.main()
