"""scripts/array_task.sh (one slurm array task) and submit_arrays.py's commands, without slurm.

Temporary data folders via APR_DATA only (never data/); nothing is submitted (dry runs).
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from datapoints import (COMPLEX_PDB, DONE_MARKER, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, TOPOLOGY,
                        TRAJECTORY, WINDOWS)
from test_main_guard import temporary_project
from test_run import FAKE_PMEMD, set_up_point

TOP = TOPOLOGY[:-len(".prmtop")]
# a "python" that only records how it was called (and so also proves when it was not called)
FAKE_PYTHON = '#!/bin/bash\necho "$*" >> "$FAKE_DIR/python_calls"\n'


class ArrayTaskTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.data = temporary_project(self.root)
        self.fake = self.root / "fake"
        (self.fake / "bin").mkdir(parents=True)
        for name, text in (("pmemd", FAKE_PMEMD), ("python", FAKE_PYTHON)):
            (self.fake / "bin" / name).write_text(text)
            (self.fake / "bin" / name).chmod(0o755)
        self.env = {**os.environ, "PATH": f"{self.fake / 'bin'}:{os.environ['PATH']}",
                    "FAKE_DIR": str(self.fake), "SLURM_ARRAY_TASK_ID": "0", "SLURM_CPUS_PER_TASK": "3"}

    def tearDown(self):
        self._tmp.cleanup()

    def worklist(self, *lines):
        path = self.root / "worklist.txt"
        path.write_text("\n".join(lines) + "\n")
        return str(path)

    def task(self, *args):
        result = subprocess.run(["bash", str(ROOT / "scripts" / "array_task.sh"), *args], env=self.env,
                                capture_output=True, text=True, timeout=120)
        return result.returncode, result.stdout + result.stderr

    def python_calls(self):
        path = self.fake / "python_calls"
        return path.read_text().splitlines() if path.exists() else []

    def test_a_simulation_task_runs_its_windows_without_python(self):
        set_up_point(self.data, "iso_00001")
        env_file = self.data / "inputs" / "run.env"
        env_file.write_text(f'DATA_DIR="{self.data}"\nPMEMD="pmemd"\nINPUT_DIR="{self.data / "inputs"}"\n'
                            f'TOPOLOGY="{TOP}"\nDONE_MARKER="{DONE_MARKER}"\nRANKS="1"\n')
        worklist = self.worklist("iso_00001/a000", "iso_00001/a001", "iso_00001/a002", "iso_00001/a003")
        code, out = self.task("sim", worklist, "2", "1", str(env_file))      # task 0 + offset 1: lines 3-4
        self.assertEqual(code, 0, out)
        self.assertIn("2 window(s) to run, 3 at a time (3 CPU(s)", out)
        windows = self.data / "iso_00001" / "windows"
        self.assertEqual([w for w in ("a000", "a001", "a002", "a003") if (windows / w / DONE_MARKER).exists()],
                         ["a002", "a003"])
        self.assertEqual(self.python_calls(), [])                         # no Python at all

    def test_setup_analysis_and_implicit_tasks_hand_their_slice_and_flags_to_the_scripts(self):
        names = self.worklist("iso_00001", "iso_00002", "iso_00003")
        self.task("setup", names, "2", "0")
        self.task("analysis", names, "2", "1", "--clean-after-analysis")
        self.task("implicit", names, "1", "1", "--p000", "--clean")
        calls = self.python_calls()
        self.assertEqual(calls[0], f"{ROOT / 'scripts' / 'main.py'} iso_00001 iso_00002")
        self.assertEqual(calls[1], f"{ROOT / 'scripts' / 'analysis.py'} --no-collect --clean-after-analysis iso_00003")
        self.assertEqual(calls[2], f"{ROOT / 'scripts' / 'run_implicit.py'} --cpus 3 --p000 --clean iso_00002")

    def test_a_task_past_the_end_of_the_list_does_nothing(self):
        code, out = self.task("setup", self.worklist("iso_00001"), "1", "5")
        self.assertEqual(code, 0, out)
        self.assertIn("0 item(s)", out)
        self.assertEqual(self.python_calls(), [])


class SubmitCommandTests(unittest.TestCase):
    """Dry runs of submit_arrays.py: the sbatch commands it would use."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.data = temporary_project(self.root)
        self.env = {**os.environ, "APR_DATA": str(self.data)}

    def tearDown(self):
        self._tmp.cleanup()

    def submit_arrays(self, *args):
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / "submit_arrays.py"), *args],
                                env=self.env, capture_output=True, text=True, timeout=120)
        return result.returncode, result.stdout + result.stderr

    def sbatch_lines(self, out):
        return [line.strip() for line in out.splitlines() if line.strip().startswith("sbatch ")]

    def parametrized(self, name):
        dp = self.data / name
        dp.mkdir()
        for filename in (GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD):
            (dp / filename).write_text("x\n")
        return dp

    def test_simulation_tasks_get_run_env(self):
        set_up_point(self.data, "iso_00001", finished=WINDOWS[:-3])
        code, out = self.submit_arrays("sim")
        self.assertEqual(code, 0, out)
        self.assertIn("[sim] 3 windows to simulate", out)
        (line,) = self.sbatch_lines(out)
        self.assertRegex(line, rf"job_array\.slurm sim \S+worklist_sim_\d+_\d+\.txt \d+ 0 "
                               rf"{re.escape(str(self.data / 'inputs' / 'run.env'))}$")
        self.assertFalse((self.data / "inputs" / "run.env").exists())        # dry run writes nothing

    def test_analysis_tasks_follow_the_config_and_pass_clean(self):
        set_up_point(self.data, "iso_00001", finished=WINDOWS)
        code, out = self.submit_arrays("analysis", "--clean-after-analysis")
        self.assertEqual(code, 0, out)
        analysis, collect = self.sbatch_lines(out)
        self.assertRegex(analysis, r"job_array\.slurm analysis \S+\.txt 1 0 --clean-after-analysis$")
        self.assertNotIn("enhanced", analysis)                                 # analysis.enhanced decides
        self.assertTrue(collect.endswith("job_array.slurm collect"))

    def test_implicit_command(self):
        self.parametrized("iso_00001")
        self.parametrized("iso_00002")
        done = self.data / "iso_00002" / "windows_implicit" / "implicit_p000"
        done.mkdir(parents=True)
        (done / TRAJECTORY).write_text("x")
        (done / DONE_MARKER).write_text("")
        (self.data / "raw").mkdir()
        (self.data / "raw" / "raw.smi").write_text("CC\n")

        code, out = self.submit_arrays("implicit")
        self.assertEqual(code, 2)
        self.assertIn("choose --implicit-only, --p000 or both", out)

        code, out = self.submit_arrays("implicit", "--p000", "--clean")
        self.assertEqual(code, 0, out)
        self.assertIn("[implicit] 1 data points", out)                          # iso_00002 is done, raw not parametrized
        (line,) = self.sbatch_lines(out)
        self.assertIn("--cpus-per-task=1", line)
        self.assertRegex(line, r"job_array\.slurm implicit \S+worklist_implicit_\d+_\d+\.txt 1 0 --p000 --clean$")

        code, out = self.submit_arrays("implicit", "--implicit-only", "--p000")
        self.assertIn("[implicit] 2 data points", out)

    def test_status_shows_failed_windows_failed_analyses_and_guest_only_runs(self):
        dp = set_up_point(self.data, "iso_00001", finished=WINDOWS[:-1])
        (dp / "windows" / WINDOWS[-1] / "run_failed.txt").write_text(
            "date: 2026-11-23 10:00:00\nhost: node1\nstep: production_part1\nexit code: 1\nlast line: x\n")
        done = set_up_point(self.data, "iso_00002", finished=WINDOWS)
        (done / "analysis_failed.txt").write_text("2026-11-23 10:00:00\npymbar did not converge\n")
        p000 = done / "windows_implicit" / "implicit_p000"
        p000.mkdir(parents=True)
        (p000 / TRAJECTORY).write_text("x")
        (p000 / DONE_MARKER).write_text("")
        only = done / "windows_implicit" / "implicit_only"
        only.mkdir(parents=True)
        (only / "run_failed.txt").write_text("step: production\n")

        code, out = self.submit_arrays("status")
        self.assertEqual(code, 0, out)
        self.assertIn("1  window(s) failed last time  (run_failed.txt", out)
        self.assertIn(f"iso_00001/{WINDOWS[-1]}: production_part1, exit code 1", out)
        self.assertIn("1  analysis failed last time  (analysis_failed.txt)", out)
        self.assertIn("iso_00002: pymbar did not converge", out)
        self.assertIn("1  implicit_p000 done", out)
        self.assertIn("1  implicit_only failed last time (run_failed.txt)", out)


if __name__ == "__main__":
    unittest.main()
