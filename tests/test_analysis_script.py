"""scripts/analysis.py on a temporary project (never data/). The pAPRika calculation itself
(apr_analysis.run) is replaced by a small stand-in; tests/test_apr_analysis_run.py covers run()."""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from datapoints import COMPLEX_PDB, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, TOPOLOGY, TRAJECTORY, WINDOWS
from test_main_guard import temporary_project

FAKE_APR_ANALYSIS = '''
import json
from pathlib import Path
def run(dp, enhanced=False):
    dp = Path(dp)
    if dp.name.startswith("bad"):
        raise RuntimeError("pymbar did not converge")
    (dp / "plots").mkdir(exist_ok=True)
    (dp / "plots" / "binding_free_energy.jpg").write_text("plot")
    out = {"dG_bind": -3.5, "dG_sem": 0.25, "enhanced": enhanced}
    (dp / "results.json").write_text(json.dumps(out))
    return out
'''

WRAPPER = '''
import importlib.util, runpy, sys
spec = importlib.util.spec_from_file_location("apr_analysis", sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
sys.modules["apr_analysis"] = module
script = sys.argv[2]
sys.argv = [script, *sys.argv[3:]]
runpy.run_path(script, run_name="__main__")
'''


class AnalysisScriptTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.data = temporary_project(self.root)
        self.env = {**os.environ, "APR_CONFIG": str(self.root / "config.toml")}
        self.fake = self.root / "fake_apr_analysis.py"
        self.fake.write_text(FAKE_APR_ANALYSIS)
        self.csv = self.data / "inputs" / "results.csv"

    def tearDown(self):
        self._tmp.cleanup()

    def point(self, name, trajectories=WINDOWS, results=None, topologies=WINDOWS):
        dp = self.data / name
        dp.mkdir()
        for filename in (GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD):
            (dp / filename).write_text("x\n")
        for window in WINDOWS:
            (dp / "windows" / window).mkdir(parents=True)
            if window in topologies:
                (dp / "windows" / window / TOPOLOGY).write_text("x\n")
            if window in trajectories:
                (dp / "windows" / window / TRAJECTORY).write_text("x\n")
        if results is not None:
            (dp / "results.json").write_text(json.dumps(results))
        return dp

    def analysis(self, *args):
        result = subprocess.run(
            [sys.executable, "-c", WRAPPER, str(self.fake), str(ROOT / "scripts" / "analysis.py"), *args],
            env=self.env, capture_output=True, text=True, timeout=120)
        return result.returncode, result.stdout + result.stderr

    def csv_rows(self):
        with open(self.csv) as f:
            return {row["smiles_folder"]: row for row in csv.DictReader(f)}

    def test_without_names_only_points_that_need_analysis(self):
        self.point("todo")
        self.point("analyzed", results={"dG_bind": -1.0, "dG_sem": 0.1})
        self.point("partial", trajectories=WINDOWS[:-1])
        code, out = self.analysis()
        self.assertEqual(code, 0, out)
        self.assertIn("[1/1] todo", out)
        self.assertIn("todo done: -3.50 +/- 0.25 kcal/mol", out)
        self.assertIn("analysis finished: 1 done, 0 failed, 0 skipped", out)
        self.assertEqual(json.loads((self.data / "analyzed" / "results.json").read_text())["dG_bind"], -1.0)
        self.assertEqual(set(self.csv_rows()), {"todo", "analyzed"})

    def test_named_points_are_reanalyzed_and_unready_ones_skipped(self):
        self.point("analyzed", results={"dG_bind": -1.0, "dG_sem": 0.1})
        self.point("partial", trajectories=WINDOWS[:-1])
        self.point("cleaned", trajectories=WINDOWS[:5], topologies=WINDOWS[:1],
                   results={"dG_bind": -2.0, "dG_sem": 0.2})
        code, out = self.analysis("analyzed", "partial", "cleaned", "missing")
        self.assertEqual(code, 0, out)
        self.assertEqual(json.loads((self.data / "analyzed" / "results.json").read_text())["dG_bind"], -3.5)
        self.assertIn(f"partial SKIPPED: not all windows have {TRAJECTORY} (1 missing, e.g. {WINDOWS[-1]})", out)
        self.assertIn("cleaned SKIPPED: cleaned up after analysis", out)
        self.assertIn("missing SKIPPED: no such data point folder", out)
        self.assertIn("analysis finished: 1 done, 0 failed, 3 skipped", out)
        # results.csv is rebuilt from every results.json, cleaned-up points included
        rows = self.csv_rows()
        self.assertEqual(set(rows), {"analyzed", "cleaned"})
        self.assertEqual(float(rows["cleaned"]["dG_bind_kcal_per_mol"]), -2.0)

    def test_a_failure_is_recorded_the_loop_goes_on_and_the_exit_code_is_0(self):
        bad = self.point("bad_one")
        self.point("good")
        code, out = self.analysis("bad_one", "good")
        self.assertEqual(code, 0, out)
        self.assertIn("bad_one FAILED: pymbar did not converge", out)
        self.assertIn("good done:", out)
        self.assertIn("analysis finished: 1 done, 1 failed, 0 skipped", out)
        self.assertIn("pymbar did not converge", (bad / "analysis_failed.txt").read_text())

        # a later success removes the record
        bad.rename(self.data / "fixed")
        code, out = self.analysis("fixed")
        self.assertFalse((self.data / "fixed" / "analysis_failed.txt").exists())

    def test_enhanced_is_passed_on_and_no_collect_leaves_results_csv_alone(self):
        self.point("todo")
        code, out = self.analysis("--enhanced", "--no-collect", "todo")
        self.assertEqual(code, 0, out)
        self.assertTrue(json.loads((self.data / "todo" / "results.json").read_text())["enhanced"])
        self.assertFalse(self.csv.exists())

        code, out = self.analysis("todo")
        self.assertFalse(json.loads((self.data / "todo" / "results.json").read_text())["enhanced"])
        self.assertTrue(self.csv.exists())

    def test_list_file(self):
        for name in ("iso_00001", "iso_00002", "iso_00003"):
            self.point(name)
        todo = self.root / "todo.txt"
        todo.write_text("iso_00003\n# not iso_00002\niso_00001\n")
        code, out = self.analysis("--list", str(todo))
        self.assertIn("[1/2] iso_00003", out)
        self.assertIn("[2/2] iso_00001", out)
        self.assertFalse((self.data / "iso_00002" / "results.json").exists())



class AnalysisConfigTests(unittest.TestCase):
    def test_old_plot_keys_are_ignored_and_results_csv_sits_in_the_data_folder(self):
        import apr_config
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = temporary_project(root)
            config = root / "config.toml"
            config.write_text(config.read_text().replace(
                "[analysis]", "[analysis]\nenhanced = true\nconvergence_plot = false", 1))
            cfg = apr_config.load_config(config)
        self.assertNotIn("enhanced", cfg["analysis"])
        self.assertNotIn("convergence_plot", cfg["analysis"])
        self.assertEqual(Path(cfg["analysis"]["results_csv"]), data / "inputs" / "results.csv")


if __name__ == "__main__":
    unittest.main()
