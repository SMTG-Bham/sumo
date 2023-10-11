import os
import unittest

try:
    import configparser
except ImportError:
    import ConfigParser as configparser
try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

import matplotlib.pyplot

import sumo.plotting
import sumo.plotting.dos_plotter
from sumo.plotting.dos_plotter import get_cached_colour


class GetColourTestCase(unittest.TestCase):
    def setUp(self):
        config_path = os.path.join(ilr_files("sumo.plotting"), "orbital_colours.conf")

        self.config = configparser.ConfigParser()
        self.config.read(os.path.abspath(config_path))

    def test_get_colour_cache(self):
        """Check colour caching"""
        col1, cache = get_cached_colour("Hf", "s", cache={})
        col2, cache = get_cached_colour("Zr", "d", cache=cache)
        col3, cache = get_cached_colour("Hf", "s", cache=cache)
        self.assertEqual(col1, col3)
        self.assertNotEqual(col1, col2)

        # Try rebooting with new cache
        col4, cache = tuple(get_cached_colour("Zr", "d", cache={}))
        self.assertEqual(col1, col4)

    def test_get_colour_global_cache(self):
        """Check colour caching"""
        sumo.plotting.colour_cache.clear()
        col1, _ = get_cached_colour("Hf", "s")
        col2, _ = get_cached_colour("Zr", "d")
        col3, _ = get_cached_colour("Hf", "s")
        self.assertEqual(col1, col3)
        self.assertNotEqual(col1, col2)

        # Try rebooting with new cache
        sumo.plotting.colour_cache.clear()
        col4, _ = get_cached_colour("Zr", "d")
        self.assertEqual(col1, col4)

    def test_get_colour_type_error(self):
        """Check bogus colour info is rejected"""
        with self.assertRaises(TypeError):
            get_cached_colour("Na", "p", colours=("#aabbcc"))

    def test_get_colour_config(self):
        """Check orbital colours from config file"""
        col_O_p, _ = get_cached_colour("O", "p", colours=self.config)
        col_Re_d, _ = get_cached_colour("Re", "d", colours=self.config)

        self.assertEqual(col_O_p, "#0DB14B")
        self.assertEqual(col_Re_d, "#A154A1")

    def test_get_colour_mixed(self):
        """Check new colours drawn in correct sequence"""
        sumo.plotting.colour_cache.clear()
        with matplotlib.pyplot.style.context("ggplot"):
            _, __ = get_cached_colour("O", "p", colours=self.config)
            col_Hf_s, _ = get_cached_colour("Hf", "s", colours=self.config)
            _, __ = get_cached_colour("Re", "d", colours=self.config)
            col_Zr_d, _ = get_cached_colour("Zr", "d", colours=self.config)

            prop_cyc = matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]
            default_0 = prop_cyc[0]
            default_1 = prop_cyc[1]
            self.assertEqual(col_Hf_s, default_0)
            self.assertEqual(col_Zr_d, default_1)


if __name__ == "__main__":
    unittest.main()
