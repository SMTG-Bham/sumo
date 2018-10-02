import logging
from pymatgen.core.structure import Structure

class QuestaalInit(object):
    """Structure information: Questaal init.ext file

    Usually this will be instantiated with the
    :obj:`~sumo.io.questaal.QuestaalInit.from_file()` method.

    Args:
        lattice
        site
        spec
    """

    def __init__(self, lattice, site, spec=None, tol=1e-5):
        self.lattice = lattice
        self.site = site
        self.tol = tol

    @property
    def structure():
        """Pymatgen structure object from Questaal init file"""


        species = [entry['ATOM'] for entry in self.site]


        if 'spcgrp' in self.lattice and self.lattice['spcgrp']:
            assert('A' in self.lattice)
            if 'A' not in self.lattice:
                logging.info('Lattice vector B not given, assume equal to A')
                self.lattice['B'] = self.lattice['A']
            if 'C' not in self.lattice:
                logging.info('Lattice vector C not given, assume equal to A')
                self.lattice['C'] = self.lattice['C']
            if 'ALPHA' not in self.lattice:
                logging.info('Lattice angle ALPHA not given, assume right-angle')
                self.lattice['ALPHA'] = 90
            if 'BETA' not in self.lattice:
                logging.info('Lattice angle BETA not given, assume right-angle')
                self.lattice['BETA'] = 90
            if 'GAMMA' not in self.lattice:
                try:
                    spcgrp_number = int(self.lattice['spcgrp'])
                except ValueError:
                    spcgrp_number = 0
                if (167 < spcgrp_number < 195):
                    logging.info('Lattice angle GAMMA not given, '
                                 'hexagonal space group, assume 120')
                    self.lattice['GAMMA'] = 120
                else:
                    logging.info('Lattice angle GAMMA not given, assume right-angle')
                    self.lattice['GAMMA'] = 90
            
            if self.cartesian:
                logging.info("Warning: Cartesian positions used without "
                             "explicit lattice vectors")
                
            return Structure.from_spacegroup(
                self.lattice['spacegroup'], lattice, species, coords,
                coords_are_cartesian=self.cartesian,
                tol=self.tol)

    @staticmethod
    def from_file(filename, tol=1e-5):
        """Read QuestaalInit object from init.ext file

        Args:
            filename (:obj:`str`): Path to init.ext file

        Returns:
            :obj:`~sumo.io.questaal.QuestaalInit`"""

        pass
