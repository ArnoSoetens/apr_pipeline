"""The data folder as a project: <data>/inputs/config.toml, host files, SMILES list, logs.

Everything runs on temporary data folders via APR_DATA (never data/).
"""
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

SHOW_DIRECTORIES = ("import json, sys; sys.path.insert(0, %r); from apr_config import CFG; "
                    "print(json.dumps(CFG['directories']))" % str(ROOT / "apr_pipeline"))


class ProjectFolderTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.data = temporary_project(self.root)
        self.config = self.data / "inputs" / "config.toml"

    def tearDown(self):
        self._tmp.cleanup()

    def python(self, *args, data=None):
        env = {**os.environ, "APR_DATA": str(data or self.data)}
        result = subprocess.run([sys.executable, *args], env=env, capture_output=True, text=True, timeout=300)
        return result.returncode, result.stdout, result.stderr

    def test_every_folder_follows_from_the_data_folder(self):
        code, out, err = self.python("-c", SHOW_DIRECTORIES)
        self.assertEqual(code, 0, err)
        directories = json.loads(out)
        inputs = self.data / "inputs"
        self.assertEqual(Path(directories["data"]), self.data)
        self.assertEqual(Path(directories["inputs"]), inputs)
        self.assertEqual(Path(directories["host"]), inputs)
        self.assertEqual(Path(directories["logs"]), self.data / "logs")
        self.assertEqual(Path(directories["smiles_list"]), inputs / "missing.txt")
        self.assertEqual(Path(directories["scripts"]), ROOT / "scripts")
        self.assertEqual(Path(directories["templates"]), ROOT / "scripts" / "templates")
        self.assertTrue((self.data / "logs").is_dir())        # slurm needs it before a job starts

    def test_a_missing_config_says_where_it_is_expected(self):
        empty = self.root / "empty_data"
        empty.mkdir()
        code, out, err = self.python("-c", SHOW_DIRECTORIES, data=empty)
        self.assertNotEqual(code, 0)
        self.assertIn(f"no settings found: expected {empty / 'inputs' / 'config.toml'}", err)

    def test_an_old_directories_section_is_ignored_with_a_note(self):
        self.config.write_text('[directories]\ndata = "/somewhere/else"\n\n' + self.config.read_text())
        code, out, err = self.python("-c", SHOW_DIRECTORIES)
        self.assertEqual(code, 0, err)
        self.assertIn("[directories]", err)
        self.assertIn("is ignored", err)
        self.assertEqual(Path(json.loads(out)["data"]), self.data)

    def test_data_point_folders_come_from_the_smiles_list_in_inputs(self):
        (self.data / "inputs" / "molecules.txt").write_text(
            "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)F\nOC(=O)C(F)(F)F\n")
        self.config.write_text(self.config.read_text().replace(
            'smiles_list = "missing.txt"', 'smiles_list = "molecules.txt"'))
        code, out, err = self.python(str(ROOT / "scripts" / "submit_arrays.py"), "setup")   # dry run
        self.assertEqual(code, 0, err)
        created = sorted(p.name for p in self.data.iterdir() if p.is_dir() and p.name not in ("inputs", "logs"))
        self.assertEqual(len(created), 2)
        self.assertIn("[setup] 2 data points", out)

    def test_inputs_and_logs_are_never_data_points(self):
        # inputs/ holds BCD.pdb, which would otherwise look like a guest PDB waiting for placement
        code, out, err = self.python(str(ROOT / "apr_pipeline" / "placement.py"))
        self.assertEqual(code, 0, err)
        self.assertIn("no data point in", out)
        for name in ("inputs", "logs"):
            code, out, err = self.python(str(ROOT / "apr_pipeline" / "placement.py"), name)
            self.assertEqual(code, 2)
            self.assertIn("not a data point", err)
        code, out, err = self.python(str(ROOT / "scripts" / "submit_arrays.py"), "status")
        self.assertEqual(code, 0, err)
        self.assertIn("0 data points", out)
        self.assertNotIn("not a data point", out)
        self.assertEqual(sorted(p.name for p in self.data.iterdir()), ["inputs", "logs"])


if __name__ == "__main__":
    unittest.main()
