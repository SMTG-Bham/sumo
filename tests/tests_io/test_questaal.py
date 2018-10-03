import unittest
from os.path import abspath
from pkg_resources import Requirement, resource_filename
from sumo.io.questaal import QuestaalInit

class QuestaalInitTestCase(unittest.TestCase):
    def test_init_from_python(self):
        """Check Questaal input object"""
        # spcgroup example
        lattice = {'SPCGRP': 186, 'A': 3.18409958, 'C': 5.1551,
                   'UNITS': 'A', 'ALAT': 1}
        site = [{'ATOM': 'Zn', 'X': (0.6666670, 0.3333330, 0.5000000)},
                {'ATOM': 'O', 'X': (0.6666670, 0.3333330, 0.8803100)}]

        init_sym = QuestaalInit(lattice, site)

        structure = init_sym.structure
        self.assertFalse(init_sym.cartesian)

        # unitcell example
        lattice = {'ALAT': 1, 'UNITS': 'A',
                   'PLAT': [[1.5920500, -2.7575110, 0.0000000],
                            [1.5920500, 2.7575110, 0.0000000],
                            [0.0000000, 0.0000000, 5.1551000]]}
        site = [{'ATOM': 'Zn', 'X': (0.6666670, 0.3333330, 0.5000000)},
                {'ATOM': 'Zn', 'X': (0.3333330, 0.6666670, 1.0000000)},
                {'ATOM': 'O', 'X': (0.6666670, 0.3333330, 0.8803100)},
                {'ATOM': 'O', 'X': (0.3333330, 0.6666670, 0.3803100)}]

        init_plat = QuestaalInit(lattice, site)

        structure = init_plat.structure
        self.assertFalse(init_plat.cartesian)

    def test_init_coordinate_safety(self):
        """Check illegal Questaal input caught"""
        lattice = {'SPCGRP': 186, 'A': 3.18409958, 'C': 5.1551,
                   'UNITS': 'A', 'ALAT': 1}
        site = [{'ATOM': 'Zn', 'X': (0.6666670, 0.3333330, 0.5000000)},
                {'ATOM': 'O', 'POS': (0.6666670, 0.3333330, 0.8803100)}]

        with self.assertRaises(ValueError):
            init = QuestaalInit(lattice, site)

        site = [{'ATOM': 'Zn', 'C': (0.6666670, 0.3333330, 0.5000000)},
                {'ATOM': 'O', 'C': (0.6666670, 0.3333330, 0.8803100)}]
        with self.assertRaises(NotImplementedError):
            init = QuestaalInit(lattice, site)
