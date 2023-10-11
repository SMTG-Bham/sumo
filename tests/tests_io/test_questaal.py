import os
import unittest

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

import numpy as np
from pymatgen.core.lattice import Lattice

from sumo.io.questaal import QuestaalInit, dielectric_from_file


class QuestaalOpticsTestCase(unittest.TestCase):
    def setUp(self):
        self.bse_path = os.path.join(ilr_files("tests"), "data", "SnO2", "eps_BSE.out")

    def test_optics_from_bethesalpeter(self):
        energy, real, imag = dielectric_from_file(self.bse_path)
        self.assertEqual(len(energy), 8000)
        self.assertAlmostEqual(energy[1], 3.401875234404301e-003)
        self.assertAlmostEqual(real[2][0], 1.94823351685956)
        self.assertAlmostEqual(imag[3][2], 1.084598600094848e-002)


class QuestaalInitTestCase(unittest.TestCase):
    def setUp(self):
        self.ref_lat = np.array(
            [
                [1.59205, -2.757511, 0.0],
                [1.59205, 2.757511, 0.0],
                [0.0, 0.0, 5.1551],
            ]
        )

        self.ref_pmg_lat = Lattice(self.ref_lat)

    def test_init_from_python(self):
        """Check Questaal input object"""
        # spcgroup example
        lattice = {
            "SPCGRP": 186,
            "A": 3.18409958,
            "C": 5.1551,
            "UNITS": "A",
            "ALAT": 1,
        }
        site = [
            {"ATOM": "Zn", "X": (0.6666670, 0.3333330, 0.5000000)},
            {"ATOM": "O", "X": (0.6666670, 0.3333330, 0.8803100)},
        ]

        init_sym = QuestaalInit(lattice, site)

        init_sym.structure
        self.assertFalse(init_sym.cartesian)

        # Check ALAT ok
        init_sym_alat = QuestaalInit(
            {
                "SPCGRP": 186,
                "A": 0.318409958,
                "C": 0.51551,
                "UNITS": "A",
                "ALAT": 10,
            },
            site,
        )
        self.assertEqual(init_sym.structure, init_sym_alat.structure)

        # unitcell example
        lattice = {
            "ALAT": 1,
            "UNITS": "A",
            "PLAT": [
                [1.5920500, -2.7575110, 0.0000000],
                [1.5920500, 2.7575110, 0.0000000],
                [0.0000000, 0.0000000, 5.1551000],
            ],
        }
        site = [
            {"ATOM": "Zn", "X": (0.6666670, 0.3333330, 0.5000000)},
            {"ATOM": "Zn", "X": (0.3333330, 0.6666670, 1.0000000)},
            {"ATOM": "O", "X": (0.6666670, 0.3333330, 0.8803100)},
            {"ATOM": "O", "X": (0.3333330, 0.6666670, 0.3803100)},
        ]

        init_plat = QuestaalInit(lattice, site)

        init_plat.structure
        self.assertFalse(init_plat.cartesian)

        init_plat_alat = QuestaalInit(lattice, site)
        init_plat_alat.lattice["PLAT"] = np.array(init_plat_alat.lattice["PLAT"]) * 0.1
        init_plat_alat.lattice["ALAT"] = 10
        self.assertEqual(init_plat.structure, init_plat_alat.structure)

        # Bohr, Cartesian units
        bohr_lattice = {
            "ALAT": 1,
            "UNITS": None,
            "PLAT": [
                [3.00853848, -5.21094058, 0.0],
                [3.00853848, 5.21094058, 0.0],
                [0.0, 0.0, 9.74172715],
            ],
        }
        bohr_cart_sites = [
            {"ATOM": "Zn", "POS": (3.00853848, -1.73698367, 4.87086358)},
            {"ATOM": "Zn", "POS": (3.00853848, 1.73698367, 9.74172715)},
            {"ATOM": "O", "POS": (3.00853848, -1.73698367, 8.57573983)},
            {"ATOM": "O", "POS": (3.00853848, 1.73698367, 3.70487625)},
        ]

        bohr_init_lat = QuestaalInit(bohr_lattice, bohr_cart_sites)
        # Cartesian detected
        self.assertTrue(bohr_init_lat.cartesian)
        # unit conversion: lattice
        self.assertLess(
            (
                abs(
                    np.array(self.ref_pmg_lat.abc)
                    - np.array(bohr_init_lat.structure.lattice.abc)
                )
            ).max(),
            1e-5,
        )
        # unit conversion: sites  (Use distance matrix so translation ignored)
        self.assertLess(
            (
                abs(
                    init_plat.structure.distance_matrix
                    - bohr_init_lat.structure.distance_matrix
                )
            ).max(),
            1e-5,
        )
        # Ignore_units option
        bohr_init_noconvert = QuestaalInit(
            bohr_lattice, bohr_cart_sites, ignore_units=True
        )
        self.assertAlmostEqual(bohr_init_noconvert.structure.lattice.abc[2], 9.74172715)
        self.assertAlmostEqual(
            bohr_init_noconvert.structure.distance_matrix[0, -1], 3.66441077
        )

    def test_init_coordinate_safety(self):
        """Check illegal Questaal input caught"""
        lattice = {
            "SPCGRP": 186,
            "A": 3.18409958,
            "C": 5.1551,
            "UNITS": "A",
            "ALAT": 1,
        }
        site = [
            {"ATOM": "Zn", "X": (0.6666670, 0.3333330, 0.5000000)},
            {"ATOM": "O", "POS": (0.6666670, 0.3333330, 0.8803100)},
        ]

        with self.assertRaises(ValueError):
            QuestaalInit(lattice, site)

        site = [
            {"ATOM": "Zn", "C": (0.6666670, 0.3333330, 0.5000000)},
            {"ATOM": "O", "C": (0.6666670, 0.3333330, 0.8803100)},
        ]
        with self.assertRaises(NotImplementedError):
            QuestaalInit(lattice, site)

    def test_init_from_file(self):
        zno_path = os.path.join(ilr_files("tests"), "data", "ZnO")

        init1 = QuestaalInit.from_file(
            os.path.join(zno_path, "init.zno_nosym"), preprocessor=False
        )
        init2 = QuestaalInit.from_file(
            os.path.join(zno_path, "init.zno_sym"), preprocessor=False
        )
        nosym_structure = init1.structure
        sym_structure = init2.structure

        self.assertLess(
            (
                abs(
                    np.array(self.ref_pmg_lat.abc)
                    - np.array(nosym_structure.lattice.abc)
                )
            ).max(),
            1e-5,
        )
        self.assertLess(
            (
                abs(
                    np.array(self.ref_pmg_lat.angles)
                    - np.array(nosym_structure.lattice.angles)
                )
            ).max(),
            1e-3,
        )

        self.assertLess(
            (
                abs(
                    np.array(self.ref_pmg_lat.abc) - np.array(sym_structure.lattice.abc)
                )
            ).max(),
            1e-5,
        )
        self.assertLess(
            (
                abs(
                    np.array(self.ref_pmg_lat.angles)
                    - np.array(sym_structure.lattice.angles)
                )
            ).max(),
            1e-3,
        )
