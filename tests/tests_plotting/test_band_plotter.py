import unittest

from sumo.plotting.bs_plotter import SBSPlotter


class SanitiseLabelTestCase(unittest.TestCase):
    def test_sanitise_label(self):
        for label_in, label_out in (
            ("X", "X"),
            ("X 1", "X 1"),
            ("X @", "X "),
            ("X@", "X"),
            ("X@@@", "X"),
            ("HEX", "HEX"),
            ("HE@X", "HE@X"),
            ("HE@X 1", "HE@X 1"),
            ("@", None),
            ("@X", None),
            ("@HEX", None),
        ):
            self.assertEqual(SBSPlotter._sanitise_label(label_in), label_out)

    def test_sanitise_label_group(self):
        for label_in, label_out in (
            ("X", "X"),
            ("X 1", "X 1"),
            ("X @", "X "),
            ("X@", "X"),
            ("X@@@", "X"),
            ("HEX", "HEX"),
            ("HE@X", "HE@X"),
            ("HE@X 1", "HE@X 1"),
            ("@", None),
            ("@X", None),
            ("@HEX", None),
            (r"X$\mid$Y", r"X$\mid$Y"),
            (r"X$\mid$Y$\mid$Z", r"X$\mid$Y$\mid$Z"),
            (r"@X$\mid$Y$\mid$Z", r"Y$\mid$Z"),
            (r"X$\mid$@Y$\mid$Z", r"X$\mid$Z"),
            (r"X$\mid$@Z", r"X"),
            (r"X@$\mid$Y", r"X$\mid$Y"),
            (r"X@$\mid$Y@@", r"X$\mid$Y"),
            (r"@X@$\mid$Y@", r"Y"),
            (r"X@$\mid$@Y", r"X"),
            (r"@X@$\mid$@Y", None),
        ):
            self.assertEqual(SBSPlotter._sanitise_label_group(label_in), label_out)
