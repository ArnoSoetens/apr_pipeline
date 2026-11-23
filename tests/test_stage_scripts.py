"""Which data points a stage script runs on. Runs structure.py for real (RDKit only) on a
temporary project, never on data/."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_main_guard import temporary_project

PFBA_ANION = "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F"


class StageScriptTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.data = temporary_project(self.root)
        self.env = {**os.environ, "APR_DATA": str(self.data)}

    def tearDown(self):
        self._tmp.cleanup()

    def structure(self, *names):
        result = subprocess.run([sys.executable, str(ROOT / "apr_pipeline" / "structure.py"), *names],
                                env=self.env, capture_output=True, text=True, timeout=300)
        return result.returncode, result.stdout + result.stderr

    def smiles_point(self, name, smiles=PFBA_ANION):
        dp = self.data / name
        dp.mkdir()
        (dp / f"{name}.smi").write_text(smiles + "\n")
        return dp

    def test_without_names_only_points_that_need_the_stage(self):
        todo = self.smiles_point("todo")
        built = self.smiles_point("built")
        (built / "UNL.pdb").write_text("kept\n")             # next stage: placement

        code, out = self.structure()
        self.assertEqual(code, 0, out)
        self.assertIn("[1/1] todo", out)
        self.assertIn("structure finished: 1 done, 0 failed, 0 skipped", out)
        self.assertTrue((todo / "UNL.pdb").is_file())
        self.assertEqual((built / "UNL.pdb").read_text(), "kept\n")
        self.assertEqual(json.loads((todo / "guest_input.json").read_text())["net_charge"], -1)

        code, out = self.structure()
        self.assertEqual(code, 0, out)
        self.assertIn("no data point in", out)

    def test_by_name_the_stage_is_redone_with_a_note_about_later_stages(self):
        built = self.smiles_point("built")
        (built / "UNL.pdb").write_text("old\n")
        code, out = self.structure("built")
        self.assertEqual(code, 0, out)
        self.assertIn("note: built already has output of later stages", out)
        self.assertIn("main.py built --from structure", out)
        self.assertNotEqual((built / "UNL.pdb").read_text(), "old\n")

    def test_range_and_list_file(self):
        for name in ("iso_00001", "iso_00002", "iso_00004", "iso_00005"):
            self.smiles_point(name)
        todo = self.root / "todo.txt"
        todo.write_text("# just these two\niso_00004\n\niso_00005  # and this one\n")

        code, out = self.structure("iso_00001-iso_00003")
        self.assertEqual(code, 0, out)
        self.assertIn("no folder for 1 of 3 (iso_00003)", out)
        self.assertIn("structure finished: 2 done, 0 failed, 0 skipped", out)

        code, out = self.structure("--list", str(todo))
        self.assertEqual(code, 0, out)
        self.assertIn("[1/2] iso_00004", out)
        self.assertIn("[2/2] iso_00005", out)
        for name in ("iso_00001", "iso_00002", "iso_00004", "iso_00005"):
            self.assertTrue((self.data / name / "UNL.pdb").is_file())

        code, out = self.structure("iso_00001/a000")
        self.assertEqual(code, 2)
        self.assertIn("whole data points", out)

    def test_a_failure_is_recorded_and_gives_exit_code_1(self):
        bad = self.smiles_point("bad", "C((")
        self.smiles_point("good")
        code, out = self.structure("bad", "good")
        self.assertEqual(code, 1)
        self.assertIn("bad FAILED in structure:", out)
        self.assertIn("1 done, 1 failed", out)
        self.assertTrue((bad / "setup_failed.txt").read_text().splitlines()[1].startswith("structure: "))


if __name__ == "__main__":
    unittest.main()
