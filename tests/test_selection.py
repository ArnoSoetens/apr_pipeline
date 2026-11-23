import argparse
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apr_pipeline"))
from selection import add_selection_args, expand_range, selected


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.data = Path(self._tmp.name) / "data"
        self.data.mkdir()
        for name in ("iso_00001", "iso_00002", "iso_00003", "iso_00005", "PFBA_m_down"):
            (self.data / name).mkdir()

    def tearDown(self):
        self._tmp.cleanup()

    def select(self, argv, windows=False, todo=lambda: ["TODO"]):
        parser = argparse.ArgumentParser()
        add_selection_args(parser, windows=windows)
        args = parser.parse_args(argv)
        out = io.StringIO()
        with redirect_stdout(out):
            result = selected(args, todo=todo, data_dir=self.data, windows=windows)
        self.printed = out.getvalue()
        return result

    def list_file(self, text):
        path = Path(self._tmp.name) / "todo.txt"
        path.write_text(text)
        return str(path)

    def test_nothing_given_uses_the_scripts_default(self):
        self.assertEqual(self.select([]), ["TODO"])

    def test_plain_names_in_order_without_duplicates(self):
        self.assertEqual(self.select(["iso_00002", "PFBA_m_down", "iso_00002"]),
                         ["iso_00002", "PFBA_m_down"])

    def test_a_name_without_folder_is_kept_for_the_script_to_report(self):
        self.assertEqual(self.select(["iso_09999"]), ["iso_09999"])

    def test_range_gives_existing_folders_and_notes_the_rest(self):
        self.assertEqual(self.select(["iso_00001-iso_00005"]),
                         ["iso_00001", "iso_00002", "iso_00003", "iso_00005"])
        self.assertIn("no folder for 1 of 5 (iso_00004)", self.printed)

    def test_bad_ranges_are_refused(self):
        for entry in ("iso_00005-iso_00001", "iso_001-iso_00010"):
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                self.select([entry])
        self.assertIsNone(expand_range("iso_00001-pfas_00003"))     # not a range: a (missing) name
        self.assertEqual(expand_range("p1-p3"), ["p1", "p2", "p3"])

    def test_a_folder_that_looks_like_a_range_is_a_name(self):
        (self.data / "run1-run2").mkdir()
        self.assertEqual(self.select(["run1-run2"]), ["run1-run2"])

    def test_list_file_with_comments_blank_lines_and_ranges(self):
        path = self.list_file("# todo for monday\n\niso_00003   # redo\nPFBA_m_down iso_00001-iso_00002\n")
        self.assertEqual(self.select(["iso_00005", "--list", path]),
                         ["iso_00005", "iso_00003", "PFBA_m_down", "iso_00001", "iso_00002"])

    def test_an_empty_list_file_selects_nothing(self):
        self.assertEqual(self.select(["--list", self.list_file("# nothing yet\n")]), [])

    def test_missing_list_file(self):
        with self.assertRaises(FileNotFoundError):
            self.select(["--list", str(Path(self._tmp.name) / "nope.txt")])

    def test_windows_only_where_the_script_takes_them(self):
        with self.assertRaisesRegex(ValueError, "whole data points"):
            self.select(["iso_00001/a000"])
        self.assertEqual(self.select(["iso_00001/a000-a002", "iso_00002"], windows=True),
                         [("iso_00001", "a000"), ("iso_00001", "a001"), ("iso_00001", "a002"),
                          ("iso_00002", None)])

    def test_names_that_leave_the_data_folder_are_refused(self):
        for entry in ("..", "iso_00001/../x", "/a000"):
            with self.subTest(entry=entry), self.assertRaises(ValueError):
                self.select([entry], windows=True)


if __name__ == "__main__":
    unittest.main()
