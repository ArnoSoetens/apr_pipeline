import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import simulate_guest_batch


class GuestBatchTests(unittest.TestCase):
    def make_data_point(self, data_dir, name, complete=True):
        data_point = data_dir / name
        data_point.mkdir()
        for suffix in ("pdb", "mol2", "frcmod"):
            if complete or suffix != "frcmod":
                (data_point / f"UNL.{suffix}").write_text("input\n", encoding="utf-8")
        return data_point

    def test_batch_continues_after_failure_and_skips_incomplete_data_points(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir)
            failed = self.make_data_point(data_dir, "PFBA_m_down")
            succeeded = self.make_data_point(data_dir, "PFBA_m_up")
            skipped = self.make_data_point(data_dir, "PFBS_mt_down", complete=False)
            fake_simulator = data_dir / "fake_simulator.py"
            fake_simulator.write_text(
                "import argparse\n"
                "import sys\n"
                "from pathlib import Path\n"
                "parser = argparse.ArgumentParser()\n"
                "parser.add_argument('pdb')\n"
                "parser.add_argument('mol2')\n"
                "parser.add_argument('--frcmod')\n"
                "parser.add_argument('--output-dir')\n"
                "parser.add_argument('--pmemd-command')\n"
                "args = parser.parse_args()\n"
                "out = Path(args.output_dir)\n"
                "print('invoked')\n"
                "if Path(args.pdb).parent.name == 'PFBA_m_down':\n"
                "    (out / 'partial.out').write_text('partial')\n"
                "    sys.exit(1)\n"
                "(out / 'guest.prmtop').write_text('topology')\n"
                "(out / 'production.nc').write_text('trajectory')\n",
                encoding="utf-8",
            )

            result = simulate_guest_batch.run_batch(
                data_dir, cpus=2, simulator=fake_simulator, clean=True
            )

            self.assertEqual(result, 1)
            self.assertIn("invoked", (failed / "guest_implicit" / "batch.log").read_text())
            self.assertTrue((failed / "guest_implicit" / "partial.out").is_file())
            self.assertFalse((succeeded / "guest_implicit" / "batch.log").exists())
            self.assertEqual(
                {path.name for path in (succeeded / "guest_implicit").iterdir()},
                {"guest.prmtop", "production.nc"},
            )
            self.assertFalse((skipped / "guest_implicit").exists())

    def test_slurm_cpu_count_is_used_by_default(self):
        with patch.dict("os.environ", {"SLURM_CPUS_PER_TASK": "12"}):
            self.assertEqual(simulate_guest_batch._get_cpu_count(None), 12)

    def test_apr_p000_batch_mode_uses_separate_output_directory(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir)
            data_point = self.make_data_point(data_dir, "PFBA_m_down")
            fake_simulator = data_dir / "fake_simulator.py"
            fake_simulator.write_text("print('completed')\n", encoding="utf-8")

            result = simulate_guest_batch.run_batch(
                data_dir, cpus=1, simulator=fake_simulator, simulation="apr-p000"
            )

            self.assertEqual(result, 0)
            self.assertTrue(
                (data_point / "guest_apr_p000" / "batch.log").is_file()
            )
            self.assertFalse((data_point / "guest_implicit").exists())

    def test_clean_keeps_only_topology_and_production_trajectory_for_both_modes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir)
            guest = self.make_data_point(data_dir, "PFBA_m_down")
            apr = self.make_data_point(data_dir, "PFBS_mt_down")
            fake_simulator = data_dir / "fake_simulator.py"
            fake_simulator.write_text(
                "import argparse\n"
                "from pathlib import Path\n"
                "parser = argparse.ArgumentParser()\n"
                "parser.add_argument('pdb')\n"
                "parser.add_argument('mol2')\n"
                "parser.add_argument('--frcmod')\n"
                "parser.add_argument('--output-dir')\n"
                "parser.add_argument('--pmemd-command')\n"
                "args = parser.parse_args()\n"
                "out = Path(args.output_dir)\n"
                "topology = 'guest.prmtop' if out.name == 'guest_implicit' "
                "else 'guest_apr_p000.prmtop'\n"
                "(out / topology).write_text('topology')\n"
                "(out / 'production.nc').write_text('trajectory')\n"
                "(out / 'production.rst7').write_text('restart')\n"
                "(out / 'nested').mkdir()\n"
                "(out / 'nested' / 'generated.tmp').write_text('temporary')\n"
                "print('completed')\n",
                encoding="utf-8",
            )

            for data_point, simulation, output_name, topology_name in (
                (guest, "guest", "guest_implicit", "guest.prmtop"),
                (apr, "apr-p000", "guest_apr_p000", "guest_apr_p000.prmtop"),
            ):
                with self.subTest(simulation=simulation):
                    result = simulate_guest_batch.run_batch(
                        data_dir,
                        cpus=1,
                        simulator=fake_simulator,
                        simulation=simulation,
                        clean=True,
                    )

                    self.assertEqual(result, 0)
                    output_dir = data_point / output_name
                    self.assertEqual(
                        {path.name for path in output_dir.iterdir()},
                        {topology_name, "production.nc"},
                    )
                    self.assertTrue(
                        all((data_point / f"UNL.{suffix}").is_file()
                            for suffix in ("pdb", "mol2", "frcmod"))
                    )


if __name__ == "__main__":
    unittest.main()
