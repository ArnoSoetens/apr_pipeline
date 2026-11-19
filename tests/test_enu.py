import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "other"))

import enu


class EnuTests(unittest.TestCase):
    def test_runs_enumerator_with_defaults_and_writes_smiles(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)

            def create_enumerator_output(command, cwd, check):
                (Path(cwd) / "conf.xml").write_text(
                    "<isomers><constitutional_SMILES>CC</constitutional_SMILES>"
                    "<constitutional_SMILES>CCC</constitutional_SMILES></isomers>"
                )

            with (
                patch.object(
                    sys,
                    "argv",
                    ["enu.py", "C[1-9]H[2-20]", "--name", "conf.txt",
                     "--output_directory", str(output_dir)],
                ),
                patch.object(enu.subprocess, "run", side_effect=create_enumerator_output) as run,
            ):
                enu.main()

            self.assertEqual(
                run.call_args.args[0],
                [
                    str(enu.ENUMERATOR),
                    "-formula",
                    "C[1-9]H[2-20]",
                    "-output",
                    "conf.xml",
                    "-output_directory",
                    str(output_dir),
                    "-unsaturations",
                    "0",
                ],
            )
            self.assertEqual(run.call_args.kwargs, {"cwd": output_dir, "check": True})
            self.assertEqual((output_dir / "conf.txt").read_text(), "CC\nCCC\n")
            self.assertFalse((output_dir / "conf.xml").exists())

    def test_keep_flag_preserves_xml_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)

            def create_enumerator_output(command, cwd, check):
                (Path(cwd) / "conf.xml").write_text(
                    "<isomers><constitutional_SMILES>CC</constitutional_SMILES></isomers>"
                )

            with (
                patch.object(
                    sys,
                    "argv",
                    [
                        "enu.py",
                        "C[1-9]H[2-20]",
                        "--name",
                        "conf.txt",
                        "--output_directory",
                        str(output_dir),
                        "--keep",
                    ],
                ),
                patch.object(enu.subprocess, "run", side_effect=create_enumerator_output),
            ):
                enu.main()

            self.assertTrue((output_dir / "conf.xml").exists())


if __name__ == "__main__":
    unittest.main()
