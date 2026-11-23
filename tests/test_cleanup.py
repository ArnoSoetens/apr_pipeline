import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import cleanup
import cleanup_batch_logs
import submit_arrays
from apr_config import CFG
from datapoints import TOPOLOGY, WINDOWS, scan


class CleanupTests(unittest.TestCase):
    def test_batch_log_cleanup_removes_only_that_batches_sim_and_analysis_logs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            log_dir = root / "logs"
            log_dir.mkdir()
            data_root = root / "data"
            point = data_root / "iso_00001"
            retained_topology = point / "windows" / WINDOWS[0] / TOPOLOGY
            retained_topology.parent.mkdir(parents=True)
            retained_topology.write_text("topology")
            (point / "results.json").write_text("{}")
            targets = [
                log_dir / "sim_b02_12345_0.out",
                log_dir / "sim_b02_12345_1.out",
                log_dir / "sim_b02_23456_0.out",
                log_dir / "analysis_b02_34567_0.out",
            ]
            preserved = [
                log_dir / "sim_b01_12345_0.out",
                log_dir / "analysis_b02_34568_0.out",
                log_dir / "cleanup_logs_b02_45678.out",
                log_dir / "worklist_sim_20261027_batch_02.txt",
            ]
            for path in targets + preserved:
                path.write_text("log")

            with patch.object(cleanup_batch_logs, "DATA", data_root):
                cleanup_batch_logs.cleanup_batch_logs(
                    log_dir, 2, ["12345", "23456"], "34567", ["iso_00001"]
                )

            self.assertTrue(all(not path.exists() for path in targets))
            self.assertTrue(all(path.exists() for path in preserved))

    def test_batch_log_cleanup_preserves_logs_if_a_point_is_not_analyzed(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            log_dir = root / "logs"
            log_dir.mkdir()
            data_root = root / "data"
            point = data_root / "iso_00001"
            point.mkdir(parents=True)
            array_log = log_dir / "analysis_b02_34567_0.out"
            array_log.write_text("log")

            with patch.object(cleanup_batch_logs, "DATA", data_root):
                cleanup_batch_logs.cleanup_batch_logs(
                    log_dir, 2, [], "34567", ["iso_00001"]
                )

            self.assertTrue(array_log.exists())

    def test_standalone_cleanup_cleans_analyzed_points_in_alternate_directory(self):
        naming = CFG["naming"]
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir) / "data"
            data_dir.mkdir()
            data_point = data_dir / "iso_00001"
            windows = data_point / "windows"
            windows.mkdir(parents=True)
            (data_point / "results.json").write_text("{}")
            (data_point / "iso_00001.smi").write_text("C")
            transient = data_point / "sqm.out"
            transient.write_text("remove")
            retained = data_point / f"{naming['complex_prefix']}.pdb"
            retained.write_text("keep")

            with patch.object(
                sys,
                "argv",
                ["cleanup.py", "--data-dir", str(data_dir)],
            ):
                cleanup.main()

            self.assertFalse(transient.exists())
            self.assertTrue(retained.exists())
            self.assertTrue((data_point / "results.json").exists())

    def test_standalone_cleanup_preserves_points_without_analysis_results(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir) / "data"
            data_dir.mkdir()
            data_point = data_dir / "iso_00001"
            data_point.mkdir()
            needed = data_point / "simulation.rst7"
            needed.write_text("keep")

            with patch.object(
                sys,
                "argv",
                ["cleanup.py", "--data-dir", str(data_dir)],
            ):
                cleanup.main()

            self.assertTrue(needed.exists())

    def test_batch_split_is_balanced_and_never_empty(self):
        names = [f"iso_{index:05d}" for index in range(8)]
        batches = submit_arrays.split_batches(names, 3)

        self.assertEqual([len(batch) for batch in batches], [3, 3, 2])
        self.assertEqual([name for batch in batches for name in batch], names)
        self.assertEqual(submit_arrays.split_batches(names, 1), [names])
        self.assertEqual(submit_arrays.split_batches(names, 20), [[name] for name in names])

    def test_batch_log_cleanup_waits_for_sim_and_analysis_success(self):
        with patch.object(submit_arrays, "sbatch", return_value=["sbatch"]) as sbatch, \
                patch.object(submit_arrays, "run_sbatch", return_value="<cleanup>"):
            submit_arrays.submit_batch_log_cleanup(
                2, ["iso_00001"], ["12345", "23456"], "34567", False
            )

        self.assertEqual(sbatch.call_args.args[0], "cleanup_logs_b02")
        self.assertEqual(
            sbatch.call_args.kwargs,
            {
                "dependency": ["12345", "23456", "34567"],
                "dependency_type": "afterok",
            },
        )
        self.assertEqual(
            sbatch.call_args.args[4],
            ["cleanup_logs", 2, "12345,23456", "34567", "iso_00001"],
        )

    def test_cleaned_data_point_remains_analyzed(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_point = Path(temp_dir) / "iso_00001"
            retained_topology = data_point / "windows" / WINDOWS[0] / TOPOLOGY
            retained_topology.parent.mkdir(parents=True)
            retained_topology.write_text("topology")
            (data_point / "results.json").write_text("{}")

            self.assertEqual(scan(data_point), ("analyzed", []))

    def test_cleanup_keeps_requested_outputs_and_five_trajectories(self):
        naming = CFG["naming"]
        with tempfile.TemporaryDirectory() as temp_dir:
            data_root = Path(temp_dir)
            data_point = data_root / "iso_00001"
            windows = data_point / "windows"
            data_point.mkdir()
            windows.mkdir()
            for window in WINDOWS:
                (windows / window).mkdir()
            retained_root = [
                data_point / f"{naming['guest_prefix']}.pdb",
                data_point / f"{naming['guest_prefix']}.mol2",
                data_point / f"{naming['guest_prefix']}.frcmod",
                data_point / f"{naming['complex_prefix']}.pdb",
                data_point / "iso_00001.smi",
                data_point / "results.json",
                data_point / "binding_free_energy.jpg",
                data_point / "pull_window_histograms.jpg",
                data_point / "attach_windows_grid.jpg",
                data_point / "convergence_plot.jpg",
                data_point / "plots" / "binding_free_energy.jpg",
                data_point / "plots" / "convergence_plot.jpg",
                data_point / "analysis_failed.txt",
            ]
            retained_root += [                        # not generated by the pipeline's later steps: kept
                data_point / "guest_input.json",
                data_point / "pfba_m_down.pdb",
                data_point / f"{naming['guest_prefix']}_clash.pdb",
                data_point / "setup_failed.txt",
                data_point / "notes.txt",
                data_point / "logs" / "simulation.log",
            ]
            for path in retained_root:
                path.parent.mkdir(exist_ok=True)
                path.write_text("keep")
            setup_leftovers = [data_point / name for name in ("sqm.in", "sqm.out", "sqm.pdb", "ATOMTYPE.INF",
                                                              "ANTECHAMBER_AC.AC", "leap.log")]
            for path in setup_leftovers:
                path.write_text("remove")
            restraints = windows / "restraints.json"
            restraints.write_text("{}")

            last_attach_window = next(
                window for window in reversed(WINDOWS) if window.startswith("a")
            )
            retained_window_names = set(WINDOWS[:5]) | {last_attach_window}
            retained_trajectories = [
                windows / window / CFG["analysis"]["trajectory"]
                for window in retained_window_names
            ]
            for path in retained_trajectories:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("keep")
            removed_window = next(
                window for window in WINDOWS if window not in retained_window_names
            )
            removed_trajectory = (
                windows / removed_window / CFG["analysis"]["trajectory"]
            )
            removed_trajectory.parent.mkdir(exist_ok=True)
            removed_trajectory.write_text("remove")
            topology_name = f"{naming['complex_prefix']}-dum.prmtop"
            retained_topology = windows / "a000" / topology_name
            retained_topology.parent.mkdir(exist_ok=True)
            retained_topology.write_text("keep")
            removed_last_attach_topology = (
                windows / last_attach_window / topology_name
            )
            removed_last_attach_topology.write_text("remove")
            removed_topology = windows / "a001" / topology_name
            removed_topology.write_text("remove")
            extra_trajectory = windows / "a000" / "other.nc"
            extra_trajectory.write_text("remove")
            intermediate = data_point / "complex" / "vac.prmtop"
            intermediate.parent.mkdir()
            intermediate.write_text("remove")

            with patch.dict(CFG["directories"], data=str(data_root)):
                cleanup.cleanup_data_point(data_point)

            self.assertTrue(all(path.exists() for path in retained_root))
            self.assertTrue(restraints.exists())
            self.assertTrue(all(path.exists() for path in retained_trajectories))
            self.assertTrue(retained_topology.exists())
            self.assertFalse(any(path.exists() for path in setup_leftovers))
            self.assertFalse(removed_trajectory.exists())
            self.assertFalse(extra_trajectory.exists())
            self.assertFalse(removed_topology.exists())
            self.assertFalse(removed_last_attach_topology.exists())
            self.assertFalse(intermediate.exists())
            self.assertTrue(removed_trajectory.parent.is_dir())
            self.assertEqual(
                {path.name for path in windows.iterdir() if path.is_dir()},
                set(WINDOWS),
            )
            self.assertEqual(
                set(windows.glob(f"*/{topology_name}")),
                {retained_topology},
            )
            self.assertEqual(
                set(windows.glob(f"*/{CFG['analysis']['trajectory']}")),
                set(retained_trajectories),
            )
            self.assertFalse((data_point / "complex").exists())
            self.assertEqual(scan(data_point), ("analyzed", []))

    def analyzed_point(self, data_dir, name="iso_00001"):
        point = data_dir / name
        for window in WINDOWS:
            (point / "windows" / window).mkdir(parents=True)
            (point / "windows" / window / TOPOLOGY).write_text("t")
            (point / "windows" / window / CFG["analysis"]["trajectory"]).write_text("trajectory")
            (point / "windows" / window / "production.out").write_text("out")
        (point / "results.json").write_text("{}")
        (point / "complex").mkdir()
        (point / "complex" / "vac.prmtop").write_text("x")
        (point / "sqm.out").write_text("x")
        return point

    def snapshot(self, point):
        return sorted(str(path.relative_to(point)) for path in point.rglob("*"))

    def run_main(self, data_dir, *args):
        import io
        from contextlib import redirect_stdout
        out = io.StringIO()
        with redirect_stdout(out):
            cleanup.main(["--data-dir", str(data_dir), *args])
        return out.getvalue()

    def test_dry_run_removes_nothing_and_counts_what_a_real_run_removes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir)
            point = self.analyzed_point(data_dir)
            before = self.snapshot(point)
            plan = cleanup.cleanup_plan(point, data_dir)
            out = self.run_main(data_dir, "--dry-run", "-v")
            self.assertEqual(self.snapshot(point), before)
            self.assertIn(f"cleanup iso_00001: would remove {len(plan)} file(s)", out)
            self.assertIn("    sqm.out", out)

            out = self.run_main(data_dir)
            self.assertIn(f"cleanup iso_00001: removed {len(plan)} file(s)", out)
            self.assertTrue(all(not path.exists() for path in plan))

    def test_unanalyzed_points_and_folders_outside_the_data_dir_are_left_alone(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir) / "data"
            data_dir.mkdir()
            point = self.analyzed_point(data_dir, "iso_00002")
            (point / "results.json").unlink()
            before = self.snapshot(point)
            out = self.run_main(data_dir, "iso_00002")
            self.assertIn("SKIPPED: not analyzed", out)
            self.assertEqual(self.snapshot(point), before)

            outside = self.analyzed_point(Path(temp_dir), "elsewhere")
            with self.assertRaisesRegex(ValueError, "outside"):
                cleanup.cleanup_plan(outside, data_dir)

    def test_selection_range_and_list(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            data_dir = Path(temp_dir) / "data"
            data_dir.mkdir()
            points = [self.analyzed_point(data_dir, f"iso_0000{i}") for i in (1, 2, 3)]
            todo = Path(temp_dir) / "todo.txt"
            todo.write_text("iso_00003\n")
            self.run_main(data_dir, "iso_00001-iso_00001", "--list", str(todo))
            self.assertFalse((points[0] / "sqm.out").exists())
            self.assertTrue((points[1] / "sqm.out").exists())
            self.assertFalse((points[2] / "sqm.out").exists())


if __name__ == "__main__":
    unittest.main()
