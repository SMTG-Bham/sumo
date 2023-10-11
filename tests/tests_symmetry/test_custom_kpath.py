import os
import unittest
import warnings

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

from pymatgen.core.structure import Structure

from sumo.symmetry.custom_kpath import CustomKpath


class CustomKpathTestCase(unittest.TestCase):
    def setUp(self):
        poscar = os.path.join(ilr_files("tests"), "data", "Ge", "POSCAR")
        with warnings.catch_warnings():  # Not interested in Pymatgen warnings
            warnings.simplefilter("ignore")
            self.structure = Structure.from_file(poscar)

    def test_auto_labels(self):
        """Check default k-point labelling"""
        kpts_flat = [
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]
        ]
        kpts_multisegment = [
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.5, 0.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
        ]

        self.assertEqual(
            CustomKpath._auto_kpath_labels(kpts_flat),
            [["(1)", "(2)", "(3)", "(4)"]],
        )
        self.assertEqual(
            CustomKpath._auto_kpath_labels(kpts_multisegment),
            [["(1)", "(2)", "(3)"], ["(4)", "(5)"], ["(4)", "(1)"]],
        )

    def test_custom_path_safety(self):
        """Check incompatible paths will be caught in CustomKpath"""
        kpt_list = [
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.5, 0.0]],
            [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
        ]
        path_labels_deep_error = [["A", "B"], ["C", "D"], ["E", "F"]]
        path_labels_shallow_error = [
            ["A", "B", "C"],
            ["D", "E"],
            ["F", "G"],
            ["H", "I"],
        ]
        self.assertRaises(
            ValueError,
            CustomKpath,
            self.structure,
            kpt_list,
            path_labels=path_labels_deep_error,
        )
        self.assertRaises(
            ValueError,
            CustomKpath,
            self.structure,
            kpt_list,
            path_labels=path_labels_shallow_error,
        )

    def test_custom_path_init(self):
        """Setup a custom Kpath"""
        path_labels = [[r"\Gamma", "X", "K"], ["U", "W", "L", "K"]]
        kpt_list = [
            [[0.0, 0.0, 0.0], [0.5, 0.0, 0.5], [0.375, 0.375, 0.75]],
            [
                [0.625, 0.25, 0.625],
                [0.5, 0.25, 0.75],
                [0.5, 0.5, 0.5],
                [0.375, 0.375, 0.75],
            ],
        ]
        kpath = CustomKpath(self.structure, kpt_list, path_labels=path_labels)

        self.assertEqual(path_labels, kpath._kpath["path"])
        for label, point in {
            r"\Gamma": [0.0, 0.0, 0.0],
            "X": [0.5, 0.0, 0.5],
            "K": [0.375, 0.375, 0.75],
            "U": [0.625, 0.25, 0.625],
            "W": [0.5, 0.25, 0.75],
            "L": [0.5, 0.5, 0.5],
        }.items():
            self.assertIn(label, kpath._kpath["kpoints"])
            self.assertEqual(point, kpath._kpath["kpoints"][label])
