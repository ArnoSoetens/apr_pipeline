"""Runs scripts/main.py and the stage scripts for real, on a temporary project (never on data/)."""
import importlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "apr_pipeline"))
from datapoints import (COMPLEX_PDB, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, STAGES, TRAJECTORY,
                        HasSimulationOutput)


def temporary_project(root):
    """A data folder root/data with inputs/config.toml (a copy of the real one) and the host files.

    Use it with the environment variable APR_DATA=<the returned folder>.
    """
    data = root / "data"
    inputs = data / "inputs"
    inputs.mkdir(parents=True)
    for suffix in (".pdb", ".mol2", ".frcmod"):
        shutil.copy(ROOT / "data" / "inputs" / f"BCD{suffix}", inputs)
    config = (ROOT / "data" / "inputs" / "config.toml").read_text()
    config = re.sub(r'(?m)^smiles_list\s*=.*$', 'smiles_list = "missing.txt"', config, count=1)
    (inputs / "config.toml").write_text(config)
    return data


def snapshot(dp):
    return sorted((str(p.relative_to(dp)), p.stat().st_mtime_ns) for p in dp.rglob("*"))


class MainGuardTests(unittest.TestCase):
    def test_main_never_touches_points_with_simulation_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = temporary_project(root)
            # windows not complete, so the old main.py would have redone setup here
            dp = data / "iso_00001"
            dp.mkdir()
            (dp / "iso_00001.smi").write_text("OC(=O)C(F)(F)C(F)(F)C(F)(F)F\n")
            for name in (GUEST_PDB, COMPLEX_PDB, GUEST_MOL2, GUEST_FRCMOD):
                (dp / name).write_text("x\n")
            (dp / "windows" / "a000").mkdir(parents=True)
            (dp / "windows" / "a000" / TRAJECTORY).write_text("x\n")
            before = snapshot(dp)

            env = {**os.environ, "APR_DATA": str(data)}
            for args in ([], ["--from", "structure"], ["iso_00001"], ["--from", "apr_setup", "iso_00001"]):
                result = subprocess.run(
                    [sys.executable, str(ROOT / "scripts" / "main.py"), *args],
                    env=env, capture_output=True, text=True, timeout=300,
                )
                with self.subTest(args=args):
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("has simulation output", result.stdout)
                    self.assertIn("0 done, 0 failed, 1 skipped", result.stdout)
                    self.assertEqual(snapshot(dp), before)

            # the stage scripts: not picked without names, skipped by name
            for stage in STAGES:
                for args in ([], ["iso_00001"]):
                    result = subprocess.run(
                        [sys.executable, str(ROOT / "apr_pipeline" / f"{stage}.py"), *args],
                        env=env, capture_output=True, text=True, timeout=300,
                    )
                    with self.subTest(stage=stage, args=args):
                        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                        self.assertIn("has simulation output" if args else f"no data point in",
                                      result.stdout)
                        self.assertEqual(snapshot(dp), before)

            # and run(dp) itself refuses, for code that calls a stage directly
            for stage in STAGES:
                module = importlib.import_module(stage)
                with self.subTest(stage=stage), self.assertRaises(HasSimulationOutput):
                    module.run(dp)
            self.assertEqual(snapshot(dp), before)


if __name__ == "__main__":
    unittest.main()
