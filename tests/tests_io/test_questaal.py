import unittest
from os.path import abspath
from pkg_resources import Requirement, resource_filename
from sumo.io.questaal import QuestaalInit

class QuestaalInitTestCase(unittest.TestCase):
    def test_init_from_python(self):
        # spcgroup example
        lattice = {'SPCGRP': 186, 'A': 3.18409958, 'C': 5.1551,
                   units: 'A', alat: 1}
        site = [{'ATOM': 'Zn', 'X': (0.6666670, 0.3333330, 0.5000000)},
                {'ATOM': 'O', 'X': (0.6666670, 0.3333330, 0.8803100)}]

        init = QuestaalInit(lattice, cite)
        self.assertTrue(True)

