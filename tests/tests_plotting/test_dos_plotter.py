from __future__ import division

import unittest
from os.path import abspath
from pkg_resources import Requirement, resource_filename
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

from sumo.plotting import colour_cycle, default_colours
import sumo.plotting.dos_plotter
from sumo.plotting.dos_plotter import get_colour_for_element_and_orbital


class GetColourTestCase(unittest.TestCase):
    def setUp(self):
        # Reboot the colour cycle for consistency
        sumo.plotting.dos_plotter.col_cycle = colour_cycle()

        # Open default CLI colours config
        config_path = resource_filename(Requirement.parse('sumo'),
                                        'sumo/plotting/orbital_colours.conf')

        self.config = configparser.ConfigParser()
        self.config.read(abspath(config_path))

    def test_get_colour_cache(self):
        """Check colour caching"""
        col1 = tuple(get_colour_for_element_and_orbital('Hf', 's'))
        col2 = tuple(get_colour_for_element_and_orbital('Zr', 'd'))
        col3 = tuple(get_colour_for_element_and_orbital('Hf', 's'))

        self.assertEqual(col1, col3)
        self.assertNotEqual(col1, col2)

    def test_get_colour_type_error(self):
        """Check bogus colour info is rejected"""
        with self.assertRaises(TypeError):
            get_colour_for_element_and_orbital('Na', 'p', colours=('#aabbcc'))

    def test_get_colour_config(self):
        """Check orbital colours from config file"""
        col_O_p = get_colour_for_element_and_orbital('O', 'p',
                                                     colours=self.config)
        col_Re_d = get_colour_for_element_and_orbital('Re', 'd',
                                                      colours=self.config)

        self.assertEqual(col_O_p, '#0DB14B')
        self.assertEqual(col_Re_d, '#A154A1')

    def test_get_colour_mixed(self):
        """Check new colours drawn in correct sequence"""
        col_O_p = get_colour_for_element_and_orbital('O', 'p',
                                                     colours=self.config)
        col_Hf_s = get_colour_for_element_and_orbital('Hf', 's',
                                                      colours=self.config)
        col_Re_d = get_colour_for_element_and_orbital('Re', 'd',
                                                      colours=self.config)
        col_Zr_d = get_colour_for_element_and_orbital('Zr', 'd',
                                                      colours=self.config)

        default_0 = tuple([x / 255 for x in default_colours[0]])
        default_1 = tuple([x / 255 for x in default_colours[1]])
        self.assertEqual(tuple(col_Hf_s), default_0)
        self.assertEqual(tuple(col_Zr_d), default_1)

if __name__ == '__main__':
    unittest.main()
