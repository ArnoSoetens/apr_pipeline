"""Runs scripts/main.py for real, on a temporary project (never on data/)."""
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
from datapoints import COMPLEX_PDB, GUEST_FRCMOD, GUEST_MOL2, GUEST_PDB, TRAJECTORY


def temporary_project(root):
    """config.toml copy whose data and host folders are root/data (host files copied in)."""
    data = root / "data"
    data.mkdir()
    for suffix in (".pdb", ".mol2", ".frcmod"):
        shutil.copy(ROOT / "data" / f"BCD{suffix}", data)
    config = (ROOT / "config.toml").read_text()
    directories = "\n".join([
        "[directories]",
        f'data = "{data}"',
        f'host = "{data}"',
        f'scripts = "{ROOT / "scripts"}"',
        f'pipeline = "{ROOT / "apr_pipeline"}"',
        f'smiles_list = "{root / "missing.txt"}"',
        f'templates = "{ROOT / "scripts" / "templates"}"',
        f'inputs = "{root / "inputs"}"',
        f'logs = "{root / "logs"}"',
        "",
    ])
    config = re.sub(r"\[directories\]\n.*?(?=\n# -{10,})", directories, config, count=1, flags=re.S)
    (root / "config.toml").write_text(config)
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

            env = {**os.environ, "APR_CONFIG": str(root / "config.toml")}
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


if __name__ == "__main__":
    unittest.main()
