"""scripts/run.sh (the bash worker) and scripts/run.py (the front end), with a fake pmemd.

Everything runs on temporary data folders via APR_DATA (never data/). The real pmemd is exercised
in tests/test_setup_end_to_end.py.
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

TOP = TOPOLOGY[:-len(".prmtop")]

# A stand-in for pmemd: records its calls, how many run at once, and creates its output files.
# FAIL_IN=<window name>:<input file name> makes that step fail with an error line in its output.
FAKE_PMEMD = r'''#!/bin/bash
mkdir -p "$FAKE_DIR/running"
touch "$FAKE_DIR/running/$$"
ls "$FAKE_DIR/running" | wc -l >> "$FAKE_DIR/concurrency"
echo "$(pwd)|$*" >> "$FAKE_DIR/calls"
args=("$@"); input=""; output=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        -i) input="$2"; shift 2 ;;
        -o) output="$2"; touch "$2"; shift 2 ;;
        -r|-x|-inf|-e|-v) touch "$2"; shift 2 ;;
        *) shift ;;
    esac
done
sleep 0.2
rm -f "$FAKE_DIR/running/$$"
if [ -n "${FAIL_IN:-}" ] && [ "$(basename "$(pwd)"):$(basename "$input")" = "$FAIL_IN" ]; then
    printf 'step output\n  ERROR: vlimit exceeded\n\n' > "$output"
    exit 3
fi
'''


class FakePmemd:
    def __init__(self, root):
        self.dir = root / "fake"
        (self.dir / "bin").mkdir(parents=True)
        exe = self.dir / "bin" / "pmemd"
        exe.write_text(FAKE_PMEMD)
        exe.chmod(0o755)

    def env(self, **extra):
        return {**os.environ, "PATH": f"{self.dir / 'bin'}:{os.environ['PATH']}",
                "FAKE_DIR": str(self.dir), **extra}

    def calls(self):
        path = self.dir / "calls"
        return path.read_text().splitlines() if path.exists() else []

    def max_concurrency(self):
        return max(int(n) for n in (self.dir / "concurrency").read_text().split())


def set_up_point(data, name, windows=WINDOWS, finished=()):
    """A data point whose setup is done (fake files), optionally with finished windows."""
    dp = data / name
    dp.mkdir()
    for filename in (GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD):
        (dp / filename).write_text("x\n")
    for window in windows:
        folder = dp / "windows" / window
        folder.mkdir(parents=True)
        (folder / f"{TOP}.prmtop").write_text("x\n")
        (folder / f"{TOP}.rst7").write_text("x\n")
        if window in finished:
            (folder / TRAJECTORY).write_text("x\n")
            (folder / DONE_MARKER).write_text("")
    return dp


class RunShTests(unittest.TestCase):
    """run.sh on its own: an env file and window folders."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.fake = FakePmemd(self.root)
        self.inputs = self.root / "inputs"
        self.inputs.mkdir()
        self.write_env()
        self.windows = self.root / "iso_00001" / "windows"

    def tearDown(self):
        self._tmp.cleanup()

    def write_env(self, ranks=1):
        (self.inputs / "run.env").write_text(
            f'PMEMD="pmemd"\nINPUT_DIR="{self.inputs}"\nTOPOLOGY="{TOP}"\n'
            f'DONE_MARKER="{DONE_MARKER}"\nRANKS="{ranks}"\n')

    def window(self, name, topology=True, finished=False):
        folder = self.windows / name
        folder.mkdir(parents=True)
        if topology:
            (folder / f"{TOP}.prmtop").write_text("x\n")
            (folder / f"{TOP}.rst7").write_text("x\n")
        if finished:
            (folder / TRAJECTORY).write_text("x\n")
            (folder / DONE_MARKER).write_text("")
        return folder

    def run_sh(self, *windows, cpus=None, env=None, env_file=True, **extra):
        args = ["bash", str(ROOT / "scripts" / "run.sh")]
        if env_file:
            args += ["--env", str(self.inputs / "run.env")]
        if cpus:
            args += ["--cpus", str(cpus)]
        result = subprocess.run([*args, *map(str, windows)], env=env or self.fake.env(**extra),
                                capture_output=True, text=True, timeout=120)
        return result.returncode, result.stdout + result.stderr

    def test_windows_run_two_at_a_time_and_get_done_markers(self):
        folders = [self.window(f"a00{i}") for i in range(4)]
        code, out = self.run_sh(*folders, cpus=2)
        self.assertEqual(code, 0, out)
        self.assertIn("4 window(s) to run, 2 at a time", out)
        self.assertIn("run finished: 4 done, 0 failed, 0 skipped", out)
        self.assertEqual(len(self.fake.calls()), 12)                   # three pmemd steps each
        self.assertEqual(self.fake.max_concurrency(), 2)
        for folder in folders:
            self.assertTrue((folder / DONE_MARKER).is_file())

    def test_finished_windows_and_windows_without_topology_are_skipped(self):
        done = self.window("a000", finished=True)
        empty = self.window("a001", topology=False)
        code, out = self.run_sh(done, empty, cpus=2)
        self.assertEqual(code, 0, out)
        self.assertIn("SKIP  iso_00001/a000 (already done)", out)
        self.assertIn("SKIP  iso_00001/a001 (no topology: setup not finished)", out)
        self.assertIn("run finished: 0 done, 0 failed, 2 skipped", out)
        self.assertEqual(self.fake.calls(), [])

    def test_a_failed_step_is_recorded_and_the_other_windows_go_on(self):
        bad, good = self.window("p000"), self.window("p001")
        code, out = self.run_sh(bad, good, cpus=2, FAIL_IN="p000:production_nopbc1.in")
        self.assertEqual(code, 1, out)
        self.assertIn("FAILED iso_00001/p000 in production_part1 (exit code 3", out)
        self.assertIn("run finished: 1 done, 1 failed, 0 skipped", out)
        record = (bad / "run_failed.txt").read_text()
        self.assertIn("step: production_part1", record)
        self.assertIn("exit code: 3", record)
        self.assertIn("last line of production.out:   ERROR: vlimit exceeded", record)
        self.assertFalse((bad / DONE_MARKER).exists())
        self.assertTrue((good / DONE_MARKER).exists())

        code, out = self.run_sh(bad, good, cpus=2)                    # fixed: run it again
        self.assertEqual(code, 0, out)
        self.assertIn("run finished: 1 done, 0 failed, 1 skipped", out)
        self.assertFalse((bad / "run_failed.txt").exists())

    def test_mpi_ranks_lower_the_number_of_windows_at_once(self):
        self.write_env(ranks=2)
        folders = [self.window(f"a00{i}") for i in range(4)]
        code, out = self.run_sh(*folders, cpus=4)
        self.assertEqual(code, 0, out)
        self.assertIn("2 at a time (4 CPU(s), 2 per window", out)
        self.assertEqual(self.fake.max_concurrency(), 2)

    def test_without_settings_file_it_stops(self):
        code, out = self.run_sh(self.window("a000"), env_file=False)
        self.assertEqual(code, 2)
        self.assertIn("settings file not found", out)


class RunPyTests(unittest.TestCase):
    """run.py: which windows, the files in data/inputs/, and the hand-over to run.sh."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.data = temporary_project(self.root)
        self.fake = FakePmemd(self.root)
        self.env = self.fake.env(APR_DATA=str(self.data))

    def tearDown(self):
        self._tmp.cleanup()

    def run_py(self, *args):
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / "run.py"), *args],
                                env=self.env, capture_output=True, text=True, timeout=300)
        return result.returncode, result.stdout + result.stderr

    def ran(self):
        """(point, window) pairs that got their first pmemd step."""
        pairs = []
        for line in self.fake.calls():
            folder, arguments = line.split("|", 1)
            if "minimize_nopbc.in" in arguments:
                pairs.append((Path(folder).parents[1].name, Path(folder).name))
        return sorted(pairs)

    def test_default_runs_every_unfinished_window_of_set_up_points(self):
        set_up_point(self.data, "iso_00001", finished=WINDOWS[:-2])     # two windows left
        set_up_point(self.data, "iso_00002", finished=WINDOWS)          # all done
        (self.data / "iso_00003").mkdir()                                # not set up
        (self.data / "iso_00003" / "iso_00003.smi").write_text("CC\n")
        code, out = self.run_py("--cpus", "2")
        self.assertEqual(code, 0, out)
        self.assertEqual(self.ran(), [("iso_00001", w) for w in sorted(WINDOWS[-2:])])
        self.assertIn("2 window(s) in 1 data point(s)", out)

    def test_points_windows_ranges_and_list(self):
        set_up_point(self.data, "iso_00001", finished=["a001"])
        set_up_point(self.data, "iso_00002")
        todo = self.root / "todo.txt"
        todo.write_text("# just one\niso_00002/p000\n")
        code, out = self.run_py("--dry-run", "iso_00001/a000-a002", "--list", str(todo))
        self.assertEqual(code, 0, out)
        self.assertIn("iso_00001/a001 SKIPPED: already done", out)
        self.assertIn("iso_00001: 2 window(s): a000 a002", out)
        self.assertIn("iso_00002: 1 window(s): p000", out)

    def test_unknown_window_names_and_unready_points(self):
        set_up_point(self.data, "iso_00001", finished=WINDOWS)
        (self.data / "raw").mkdir()
        (self.data / "raw" / "guest.pdb").write_text("x\n")
        code, out = self.run_py("iso_00001/a999")
        self.assertNotEqual(code, 0)
        self.assertIn("not a configured window: a999", out)
        code, out = self.run_py("iso_00001", "raw", "missing", "raw/a000-a002")
        self.assertEqual(code, 0, out)
        self.assertEqual(out.count("raw SKIPPED"), 1)                   # once per point, not per window
        self.assertIn("iso_00001 SKIPPED: all windows already done", out)
        self.assertIn("raw SKIPPED: setup not finished (next stage: placement)", out)
        self.assertIn("missing SKIPPED: no such data point folder", out)
        self.assertIn("nothing to run", out)

    def test_dry_run_writes_nothing(self):
        set_up_point(self.data, "iso_00001")
        code, out = self.run_py("--dry-run")
        self.assertEqual(code, 0, out)
        self.assertIn(f"{len(WINDOWS)} window(s) in 1 data point(s)", out)
        self.assertEqual(self.fake.calls(), [])
        self.assertFalse((self.data / "inputs" / "run.env").exists())

    def test_a_run_writes_its_inputs_and_hands_the_windows_to_run_sh(self):
        dp = set_up_point(self.data, "iso_00001")
        config = self.data / "inputs" / "config.toml"
        config.write_text(re.sub(r"(?m)^max_jobs\s*=\s*\d+", "max_jobs = 3", config.read_text()))
        code, out = self.run_py("iso_00001/a000-a001")
        self.assertEqual(code, 0, out)
        self.assertIn("3 at a time (3 CPU(s)", out)                      # max_jobs is the default
        self.assertIn("run finished: 2 done, 0 failed, 0 skipped", out)
        inputs = self.data / "inputs"
        for name in ("minimize_nopbc.in", "production_nopbc1.in", "production_nopbc2.in", "run.env"):
            self.assertTrue((inputs / name).is_file(), name)
        self.assertIn('PMEMD="pmemd"', (inputs / "run.env").read_text())
        self.assertTrue((dp / "windows" / "a001" / DONE_MARKER).is_file())
        first_call = self.fake.calls()[0].split("|")[1]
        self.assertIn(f"-i {inputs / 'minimize_nopbc.in'}", first_call)

    def test_a_failed_window_gives_exit_code_1(self):
        set_up_point(self.data, "iso_00001")
        self.env["FAIL_IN"] = "a000:minimize_nopbc.in"
        code, out = self.run_py("--cpus", "1", "iso_00001/a000")
        self.assertEqual(code, 1, out)
        self.assertTrue((self.data / "iso_00001" / "windows" / "a000" / "run_failed.txt").is_file())


class AnalyzedPointsTests(unittest.TestCase):
    """A cleaned-up point keeps a few trajectories but no done markers: never rerun unless asked."""

    def test_analyzed_points_are_only_run_when_a_window_is_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = temporary_project(root)
            fake = FakePmemd(root)
            env = fake.env(APR_DATA=str(data))
            dp = set_up_point(data, "iso_00001", finished=WINDOWS)
            (dp / "results.json").write_text("{}\n")
            for window in WINDOWS:                                # what cleanup leaves behind
                folder = dp / "windows" / window
                (folder / DONE_MARKER).unlink()
                if window != WINDOWS[0]:
                    (folder / f"{TOP}.prmtop").unlink()
                    if window not in WINDOWS[:5]:
                        (folder / TRAJECTORY).unlink()

            def run_py(*args):
                result = subprocess.run([sys.executable, str(ROOT / "scripts" / "run.py"), *args], env=env,
                                        capture_output=True, text=True, timeout=120)
                return result.returncode, result.stdout + result.stderr

            code, out = run_py()
            self.assertIn("nothing to run", out)
            code, out = run_py("iso_00001")
            self.assertIn("iso_00001 SKIPPED: already analyzed; name a window (iso_00001/a000)", out)
            self.assertEqual(fake.calls(), [])

            code, out = run_py("iso_00001/a000")                  # asked for explicitly: it runs
            self.assertEqual(code, 0, out)
            self.assertIn("run finished: 1 done", out)


class MainWithWindowEntriesTests(unittest.TestCase):
    def test_window_entries_select_their_data_point(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = temporary_project(Path(tmp))
            set_up_point(data, "iso_00001")
            result = subprocess.run([sys.executable, str(ROOT / "scripts" / "main.py"), "iso_00001/a000-a002"],
                                    env={**os.environ, "APR_DATA": str(data)},
                                    capture_output=True, text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("[1/1] iso_00001", result.stdout)
            self.assertIn("setup is complete", result.stdout)


if __name__ == "__main__":
    unittest.main()
