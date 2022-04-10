import unittest
from os.path import join as path_join

import numpy as np
from phonopy import Phonopy
from pkg_resources import resource_filename
from pymatgen.io.vasp.inputs import Poscar

import sumo.phonon.phonopy


class PhonopyTestCase(unittest.TestCase):
    def setUp(self):
        self.phonon_data = resource_filename(
            __name__, path_join("..", "data", "RbSnI6", "phonopy", "FORCE_SETS")
        )
        poscar_path = resource_filename(
            __name__, path_join("..", "data", "RbSnI6", "phonopy", "POSCAR")
        )
        phonon_poscar = Poscar.from_file(poscar_path)
        self.phonon_structure = phonon_poscar.structure

    def test_phonopy_load(self):
        phonopy_obj = sumo.phonon.phonopy.load_phonopy(
            self.phonon_data, self.phonon_structure, np.diag([3, 3, 2])
        )
        self.assertIsInstance(phonopy_obj, Phonopy)

        self.assertEqual(phonopy_obj.force_constants.shape, (324, 324, 3, 3))
        self.assertAlmostEqual(
            phonopy_obj.force_constants.trace().trace(), 2654.3381334806545
        )
        self.assertAlmostEqual(
            phonopy_obj.force_constants[4, 3, 2, 1], 0.000552387416908992
        )
