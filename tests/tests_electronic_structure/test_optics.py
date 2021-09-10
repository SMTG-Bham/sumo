import json
import unittest
from os.path import join as path_join

import numpy as np
from numpy.testing import assert_almost_equal
from pkg_resources import resource_filename
from pymatgen.io.vasp import Vasprun

from sumo.electronic_structure.optics import calculate_dielectric_properties, kkr


class AbsorptionTestCase(unittest.TestCase):
    def setUp(self):
        diel_path = resource_filename(
            __name__, path_join("..", "data", "Ge", "ge_diel.json")
        )
        with open(diel_path) as f:
            self.ge_diel = json.load(f)

        absorption_path = resource_filename(
            __name__, path_join("..", "data", "Ge", "ge_alpha.json")
        )
        with open(absorption_path) as f:
            self.ge_abs = json.load(f)

    def test_absorption(self):
        energy, properties = calculate_dielectric_properties(
            self.ge_diel,
            {"absorption"},
        )
        self.assertIsNone(assert_almost_equal(properties["absorption"], self.ge_abs))


class KramersKronigTestCase(unittest.TestCase):
    def setUp(self):
        ge_vasprun_path = resource_filename(
            __name__, path_join("..", "data", "Ge", "vasprun.xml.gz")
        )
        self.ge_vasprun = Vasprun(ge_vasprun_path)

        self.ge_text_file = resource_filename(
            __name__, path_join("..", "data", "Ge", "optics.txt")
        )

    def test_kkr(self):
        energy, eps_real, eps_imag = self.ge_vasprun.dielectric

        de = (energy[10] - energy[0]) / 10

        def symmetrise(a):
            """Convert XX YY ZZ XY YZ XZ array to a symmetrical 3x3 matrix"""
            return [[a[0], a[3], a[5]], [a[4], a[1], a[4]], [a[5], a[4], a[2]]]

        eps_imag_3x3 = [symmetrise(a) for a in eps_imag]
        eps_real_3x3 = np.array([symmetrise(a) for a in eps_real])

        # Check difference between eps_real reported by Vasp and determined
        # by Kramers-Kronig transformation of eps_im.
        #
        # Some discrepancy is normal, check RMS is as expected
        # This is likely due to the limited precision available in vasprun

        error = kkr(de, eps_imag_3x3) - eps_real_3x3
        error_fracs = [
            eps / eps_ref
            for eps, eps_ref in zip(error.flatten(), eps_real_3x3.flatten())
            if eps_ref > 1e-2
        ]  # Exclude low-precision cases

        self.assertLess(np.sqrt((np.array(error_fracs) ** 2).mean()), 0.1)
