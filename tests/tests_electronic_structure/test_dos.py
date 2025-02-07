import os
import unittest

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

from sumo.electronic_structure.dos import load_dos
from pymatgen.electronic_structure.core import Spin


class DosTestCase(unittest.TestCase):
    def setUp(self):
        # SOC DOS calculation:
        self.vr_path = os.path.join(ilr_files("tests"), "data", "Cs2SnBr6", "vasprun.xml.gz")

    def test_load_dos(self):
        dos, pdos = load_dos(self.vr_path)  # previously would fail due to old SOC vr handling
        self.assertEqual(len(dos.energies), 2000)  # NEDOS
        self.assertEqual(list(dos.densities.keys()), [Spin.up])  # no Spin down key
        self.assertEqual(len(dos.densities[Spin.up]), len(dos.energies))  # matching array lengths
