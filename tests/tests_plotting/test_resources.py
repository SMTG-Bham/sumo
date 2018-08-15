import unittest
from os.path import isfile
from pkg_resources import Requirement, resource_filename
from sumo.plotting import sumo_base_style, sumo_dos_style, sumo_bs_style
from sumo.plotting import sumo_phonon_style, sumo_optics_style

class StyleSheetsTestCase(unittest.TestCase):
    def test_style_files(self):
        """Check style sheets exist"""
        for style in (sumo_base_style,
                      sumo_dos_style,
                      sumo_bs_style,
                      sumo_phonon_style,
                      sumo_optics_style):
            self.assertTrue(isfile(style))
