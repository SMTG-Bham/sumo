import logging
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

_bohr_to_angstrom = 0.5291772

class QuestaalInit(object):
    """Structure information: Questaal init.ext file

    Usually this will be instantiated with the
    :obj:`~sumo.io.questaal.QuestaalInit.from_file()` method.
    Data from each file section is given as a separate input argument and
    stored as a property.

    Args:
        lattice (:obj:`dict`):
            There are two main forms of lattice data (stored in self.lattice):
            - Explicit lattice vectors expressed with PLAT
            - Lattice angles/lengths expressed as A, B, C, ALPHA, BETA, GAMMA
              - This style of input must be used with SPCGRP, specifying
                symmetry
              - Angles are given in degrees

            Examples:

                lattice = {'SPCGRP': 186, 'A': 3.18409958, 'C': 5.1551,
                           'UNITS': 'A', 'ALAT': 1}

                lattice = {'ALAT': 1, 'UNITS': 'A',
                           'PLAT': [[1.59, -2.75, 0.],
                                    [1.59, 2.75, 0.],
                                    [0., 0., 5.16]]}
        site (:obj:`list`): Site species, coordinate type and location as a
            list of dictionaries of form
            ``{'ATOM': el, COORD_TYPE: (a, b, c)}``
            where *el* is the species label, COORD_TYPE is ``'POS'``
            (Cartesian) or ``'X'`` (direct/fractional) and (a, b, c) is a tuple
            giving the site position in Cartesian or fractional coordinates.

        spec
    """

    def __init__(self, lattice, site, spec=None, tol=1e-5):
        self.lattice = lattice
        self.site = site
        self.tol = tol

        cartesian_sites = any('POS' in item for item in site)
        fractional_sites = any('X' in item for item in site)
        c_sites = any('C' in item for item in site)

        if c_sites:
            raise NotImplementedError('C position option for Questaal input '
                                      '(conventional lattice vector fractions)'
                                      ' not implemented. Life is too short!')
        if cartesian_sites and fractional_sites:
            raise ValueError("Cannot mix direct and Cartesian input")
        else:
            self.cartesian = cartesian_sites

    @property
    def structure(self):
        """Pymatgen structure object from Questaal init file

        """

        if 'SPCGRP' in self.lattice and self.lattice['SPCGRP']:
            return self._get_structure_from_spcgrp()
        else:
            return self._get_structure_from_lattice()

    def _get_species_coords(self):
        species = [entry['ATOM'] for entry in self.site]

        if self.cartesian:
            coords = [entry['POS'] for entry in self.site]
        else:
            coords = [entry['X'] for entry in self.site]

        return species, coords

    def _get_structure_from_spcgrp(self):
        assert('A' in self.lattice)
        if 'B' not in self.lattice:
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
                spcgrp_number = int(self.lattice['SPCGRP'])
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

        if 'UNITS' not in self.lattice or self.lattice['UNITS'] is None:
            for length in ('A', 'B', 'C'):
                self.lattice[length] *= _bohr_to_angstrom

        assert('ALAT' in self.lattice)
        for length in ('A', 'B', 'C'):
            self.lattice[length] *= self.lattice['ALAT']

        lattice = Lattice.from_parameters(self.lattice['A'],
                                          self.lattice['B'],
                                          self.lattice['C'],
                                          self.lattice['ALPHA'],
                                          self.lattice['BETA'],
                                          self.lattice['GAMMA'])

        species, coords = self._get_species_coords()

        return Structure.from_spacegroup(
            self.lattice['SPCGRP'], lattice, species, coords,
            coords_are_cartesian=self.cartesian,
            tol=self.tol)

    def _get_structure_from_lattice(self):
        lattice = Lattice(self.lattice['PLAT'])
        lattice = Lattice(lattice.matrix * self.lattice['ALAT'])

        species, coords = self._get_species_coords()

        return Structure(lattice, species, coords,
                         coords_are_cartesian=self.cartesian)

    @staticmethod
    def from_file(filename, preprocessor=True, tol=1e-5):
        """Read QuestaalInit object from init.ext file

        Args:
            filename (:obj:`str`): Path to init.ext file
            preprocessor (:obj:`bool`): Process file with ``rdfile`` (must be
                available on shell PATH).
            tol (:obj:`float`, optional): tolerance for symmetry operations

        Returns:
            :obj:`~sumo.io.questaal.QuestaalInit`"""

        if preprocessor:
            from subprocess import Popen, PIPE
            process = Popen(['rdfile', filename], stdout=PIPE)
            lines = process.stdout.readlines()
            #  Need to decode from bytes. Hard-coding ASCII here - it doesn't
            #  seem likely that Questaal would support unicode?
            lines = [line.decode('ascii') for line in lines]
        else:
            with open(filename, 'r') as f:
                lines = f.readlines()

        categories = {'LATTICE', 'SITE', 'SPEC'}

        # Find which lines begin a new category
        cat_lines = []
        for i, line in enumerate(lines):
            if line.strip().split()[0] in categories:
                cat_lines.append(i)
        cat_lines.append(None)  # None allows us to slice up to the file end

        # Grab the lines corresponding to each section and collect by category
        grouped_lines = {}
        for i in range(len(cat_lines) - 1):
            category = lines[cat_lines[i]].split()[0]
            grouped_lines[category] = lines[cat_lines[i]:cat_lines[i + 1]]

        # Initial cleanup: - Remove leading/trailing whitespace
        #                  - drop lines beginning with '#'
        #                  - remove category name from first line

        for category, lines in grouped_lines.items():
            lines = [line.strip() for line in lines if line.strip()[0] != '#']

            category_line_remainder = lines[0][len(category):].strip()
            lines = [category_line_remainder] + lines[1:]

            grouped_lines[category] = lines

        # Join lines and split into tags
        init_data = {}
        for category, lines in grouped_lines.items():
            tag_text = ' '.join(lines)

            if category == 'SITE':
                pass
            else:

            # We have something like this:
            # "TAG1= X Y Z TAG2 = A TAG3 =A B C D E F"
            # Split on = sign.
            # For each split, next tag symbol is at the end
            # (e.g. " X Y Z ->TAG2<-")
            tag_symbols = [x.split()[-1] for x in tag_text.split('=')]
            # Except for the last one (A B C D E ->F<-) so drop that
            tag_symbols = tag_symbols[:-1]
            # For each split, tag content is right of split up to penultimate
            # item (e.g. ->X Y Z<- TAG2)
            tag_content = [x.split()[:-1] for x in tag_text.split('=')]
            # Except for the last one, which is whole item (->A B C D E F<-)
            tag_content[-1] = tag_text.split('=')[-1].split()

            init_data[category] = dict(zip(tag_symbols, tag_content))

        print(init_data)
