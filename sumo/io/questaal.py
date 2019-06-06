from math import ceil
from itertools import chain, product
import errno
import os.path
from os import makedirs
from subprocess import Popen, PIPE
from shutil import which
import re
import logging
import numpy as np
from monty.io import zopen
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.dos import Dos, CompleteDos

from sumo.electronic_structure.optics import kkr
import sumo.electronic_structure.dos

_bohr_to_angstrom = 0.5291772
_ry_to_ev = 13.605693009


class QuestaalSite(object):
    """Structure information: Questaal site.ext file

    Usually this will be instantiated with the
    :obj:`~sumo.io.questaal.QuestaalSite.from_file()` method.

    """

    def __init__(self, nbas, vn=3., io=15, alat=1., xpos=True, read='fast',
                 sites=None, plat=[1, 0, 0, 0, 1, 0, 0, 0, 1]):
        sites = sites or []
        if nbas != len(sites):
            raise AssertionError()
        if len(plat) != 9:
            raise AssertionError()

        if read != 'fast':
            raise Exception("Algebraic expressions not supported, use 'fast'")
        if io != 15:
            raise Exception("Only site.ext format 15 supported at present")

        self.nbas, self.vn, self.io, self.alat = nbas, vn, io, alat
        self.xpos, self.read, self.sites, self.plat = xpos, read, sites, plat

        is_empty = re.compile('E\d*$')
        empty_sites = [site for site in sites
                       if is_empty.match(site['species']) is not None]
        self.nbas_empty = len(empty_sites)

    @property
    def structure(self):
        # Get lattice vectors in Angstrom
        lattice = Lattice(self.plat)
        lattice = Lattice(lattice.matrix * self.alat * _bohr_to_angstrom)

        # Get corresponding lists of species and positions by making a list of
        # pairs and unpacking with zip
        if self.xpos:
            species_coords = [(site['species'], site['pos'])
                              for site in self.sites]
            species, coords = zip(*species_coords)

            return Structure(lattice, species, coords,
                             coords_are_cartesian=False)
        else:
            species_coords = [(site['species'],
                               [x * self.alat * _bohr_to_angstrom
                                for x in site['pos']])
                              for site in self.sites]
            species, coords = zip(*species_coords)

            return Structure(lattice, species, coords,
                             coords_are_cartesian=True)

    @classmethod
    def from_file(cls, filename):
        with zopen(filename, 'rt') as f:
            lines = f.readlines()

        header = lines[0]
        sites = [line for line in lines if line[0] not in '#%']

        # Some of the header info does not use '=' so handle separately
        header_items = header.strip().split()
        if header_items[0] != '%' or header_items[1] != 'site-data':
            raise AssertionError()

        xpos = True if 'xpos' in header_items else False
        read = 'fast' if 'fast' in header_items else False

        header_clean = ' '.join(x for x in header_items if x not in
                                ('%', 'site-data', 'xpos', 'fast'))

        tags = re.findall(r'(\w+)\s*=', header_clean)  # Find tags
        # Split on tags to find tag parameters
        tag_data = re.split(r'\s*\w+\s*=\s*', header_clean)[1:]
        tag_dict = dict(zip(tags, tag_data))

        vn = float(tag_dict['vn']) if 'vn' in tag_dict else 3.
        io = int(tag_dict['io']) if 'io' in tag_dict else 15.
        nbas = int(tag_dict['nbas']) if 'nbas' in tag_dict else 15.
        alat = float(tag_dict['alat']) if 'alat' in tag_dict else 1.
        plat = ([float(x) for x in tag_dict['plat'].split()]
                if 'plat' in tag_dict else [1, 0, 0, 0, 1, 0, 0, 0, 1])

        # Convert sites to structured format
        # (If needed, could support other 'io' options here)
        sites = [{'species': site.split()[0],
                  'pos': [float(x) for x in site.split()[1:4]]}
                 for site in sites]

        return cls(nbas, vn=vn, io=io, alat=alat, xpos=xpos, read=read,
                   sites=sites, plat=plat)


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

        spec (:obj:`dict`, optional): Content of SPEC section, expressed as
            Python dict

        tol (:obj:`float`, optional): Tolerance in units of distance used when
            constructing cell from space group

        ignore_units (:obj:`bool`, optional): If True, no unit conversions will
            be applied when converting to other formats (e.g. Pymatgen
            structure). This is needed for consistent use of Bohr units when
            interpreting band structure data.

    """

    def __init__(self, lattice, site, spec=None, tol=1e-5, ignore_units=False):
        self.lattice = lattice
        self.site = site
        self.spec = spec
        self.tol = tol
        self.ignore_units = ignore_units

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
        if 'A' not in self.lattice:
            raise AssertionError()
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
                logging.info('Lattice angle GAMMA not given, '
                             'assume right-angle')
                self.lattice['GAMMA'] = 90

        if self.cartesian:
            logging.info("Warning: Cartesian positions used without "
                         "explicit lattice vectors")

        if 'UNITS' not in self.lattice or self.lattice['UNITS'] is None:
            if self.ignore_units:
                pass
            else:
                for length in ('A', 'B', 'C'):
                    self.lattice[length] *= _bohr_to_angstrom

        if 'ALAT' not in self.lattice:
            raise AssertionError()
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

        if self.ignore_units:
            pass
        else:
            if 'UNITS' not in self.lattice or self.lattice['UNITS'] is None:
                lattice = Lattice(lattice.matrix * _bohr_to_angstrom)
                coords = [tuple([x * _bohr_to_angstrom for x in site_coords])
                          for site_coords in coords]

        return Structure(lattice, species, coords,
                         coords_are_cartesian=self.cartesian)

    def to_file(self, filename):
        """Write QuestaalInit object to init file"""
        with open(filename, 'w') as f:

            f.write('LATTICE\n')
            for key, value in self.lattice.items():
                if key == 'PLAT':
                    #  Expand nested lists to one flat list.
                    #  Yes, nested list comprehensions look weird!
                    lattice_params = [c for row in self.lattice['PLAT']
                                      for c in row]
                    # Write out as string-separated row of 9
                    f.write('    PLAT= ' + ' '.join(map(str, lattice_params)))
                    f.write('\n')
                else:
                    f.write('    {0}={1}\n'.format(key, value))

            f.write('SITE\n')
            for row in self.site:
                f.write('    ATOM={0:4s}  '.format(row['ATOM']))
                if 'POS' in self.site:
                    f.write('POS= {0:11.8f} {1:11.8f} {2:11.8f}'.format(
                        *row['POS']))
                else:
                    f.write('X= {0:11.8f} {1:11.8f} {2:11.8f}'.format(
                        *row['X']))
                for key, value in row.items():
                    if key not in ('ATOM', 'POS', 'X'):
                        f.write('  {0}= {1}'.format(key, value))
                f.write('\n')

            if self.spec is not None:
                f.write('SPEC')
                for key, value in self.spec.items():
                    f.write('  {0}= {1}\n'.format(key, value))

    @staticmethod
    def from_structure(structure):
        """Generate QuestaalInit object from pymatgen structure"""
        lattice = {'ALAT': 1, 'UNITS': 'A'}
        lattice['PLAT'] = structure.lattice.matrix

        sites = [{'ATOM': site.species_string, 'X': tuple(site.frac_coords)}
                 for site in structure.sites]
        return QuestaalInit(lattice, sites)

    @staticmethod
    def from_file(filename, preprocessor=None, tol=1e-5, ignore_units=False):
        """Read QuestaalInit object from init.ext file

        Args:
            filename (:obj:`str`): Path to init.ext file
            preprocessor (:obj:`bool`): Process file with ``rdfile`` (must be
                available on shell PATH). If None, use preprocessor where
                available. If True, an error will be raised if ``rdfile`` is
                unavailable.
            tol (:obj:`float`, optional): tolerance for symmetry operations
            ignore_units (:obj:`bool`, optional): If True, no unit conversions
                will be applied when converting to other formats (e.g. Pymatgen
                structure). This is needed for consistent use of Bohr units
                when interpreting band structure data.


        Returns:
            :obj:`~sumo.io.questaal.QuestaalInit`"""

        if preprocessor is None:
            if which('rdfile') is not None:
                preprocessor = True
            else:
                preprocessor = False

        if preprocessor:
            process = Popen(['rdfile', filename], stdout=PIPE)
            lines = process.stdout.readlines()
            #  Need to decode from bytes. Hard-coding ASCII here - it doesn't
            #  seem likely that Questaal would support unicode?
            lines = [line.decode('ascii') for line in lines]
        else:
            with zopen(filename, 'r') as f:
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
                site_data = []
                #  Split on regex: ATOM tag will be removed, species is left
                #  followed by other tags, e.g.
                #  "ATOM=Zn X = 0.0 0.0 0.5 ATOM = S  X= 0.0 0.0 0.0"
                #  is split to
                #  ['', 'Zn X = 0.0 0.0 0.5', 'S X = 0.0 0.0 0.0']
                #
                atom_entries = re.split(r'ATOM\s*=\s*', tag_text)

                for line in atom_entries[1:]:  # Drop the empty first line
                    atom = line.split()[0]
                    tag_dict = {'ATOM': atom}

                    line = line[len(atom):]    # Drop species tag from line
                    tags = re.findall(r'(\w+)\s*=', line)  # Find tags
                    # Split on tags to find tag parameters
                    tag_data = re.split(r'\s*\w+\s*=\s*', line)[1:]
                    tag_dict.update(dict(zip(tags, tag_data)))

                    # Cast coordinates to tuple
                    for key in ('POS', 'X'):
                        if key in tag_dict:
                            tag_dict[key] = tuple(map(float,
                                                      tag_dict[key].split()))

                    site_data.append(tag_dict)

                init_data['SITE'] = site_data

            else:
                float_params = ('A', 'B', 'C',
                                'ALPHA', 'BETA', 'GAMMA',
                                'ALAT')

                unsupported_params = ('GENS')

                tags = re.findall(r'(\w+)\s*=', tag_text)  # Find tags
                # Split on tags to find tag parameters
                tag_data = re.split(r'\s*\w+\s*=\s*', tag_text)[1:]
                tag_dict = dict(zip(tags, tag_data))

                if 'SPCGRP' in tag_dict:
                    try:
                        tag_dict['SPCGRP'] = int(tag_dict['SPCGRP'])
                    except ValueError:
                        pass

                if 'PLAT' in tag_dict:
                    lattice = tuple(map(float, tag_dict['PLAT'].split()))
                    if len(lattice) != 9:
                        raise AssertionError()
                    tag_dict['PLAT'] = [[lattice[0], lattice[1], lattice[2]],
                                        [lattice[3], lattice[4], lattice[5]],
                                        [lattice[6], lattice[7], lattice[8]]]

                for float_param in float_params:
                    if float_param in tag_dict:
                        tag_dict[float_param] = float(tag_dict[float_param])

                for unsupported_param in unsupported_params:
                    if unsupported_param in tag_dict:
                        raise NotImplementedError(
                            'Questaal tag {0}_{1} is not supported'.format(
                                category, unsupported_param))

                init_data[category] = tag_dict

            if 'SPEC' not in init_data or init_data['SPEC'] == {}:
                init_data['SPEC'] = None

        return QuestaalInit(init_data['LATTICE'],
                            init_data['SITE'],
                            spec=init_data['SPEC'],
                            tol=tol)


def write_kpoint_files(filename, kpoints, labels, alat=1,
                       make_folders=False, directory=None, cart_coords=False,
                       **kwargs):
    """Write syml file for Questaal kpoints

    The interface imitates the VASP KPOINTS file writer for simplicity of
    integration into Sumo, but there are some conceptual differences:

    If *labels* is None, then *kpoints* will be written to a simple file in
    "List mode".

    If labels are provided, then labelled points will be extracted from the
    list of kpoints and "line mode" is used to create a compact band-structure
    input file with the same number of k-points.

    Args:
        filename (:obj:`str`): Path to init.ext file. (Extension is used to
            name syml file).
        kpoints (:obj:`numpy.ndarray`): The k-point coordinates along the
            high-symmetry path. For example::

                [[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0], [0.5, 0, 0.25],
                [0.5, 0, 0.5]]

        labels (:obj:`list`) The high symmetry labels for each k-point (will be
            an empty :obj:`str` if the k-point has no label). For example::

                ['\Gamma', '', 'X', '', 'Y']

        make_folders (:obj:`bool`, optional): Generate folders and copy in
            required files if found from the current directory.
        directory (:obj:`str`, optional): The output file directory.
        cart_coords (:obj:`bool`, optional): Whether the k-points are returned
            in Cartesian or reciprocal coordinates. Defaults to ``False``
            (fractional coordinates). Note that Questaal uses crystallographic-
            style Cartesian coordinates and a factor of 1/(2 pi) should be
            included in the kpoints data.
    """

    for key, value in kwargs.items():
        if value is not None:
            logging.info('Ignoring k-point write option "{0}"; not '
                         'implemented for Questaal calculations.'.format(key))

    ext = filename.split('.')[-1]
    logging.info('System id from init filename: {0}'.format(ext))

    if directory is not None:
        path = directory
    else:
        path = os.path.curdir

    if make_folders:
        path = os.path.join(path, 'band-calc')

        try:
            makedirs(path)
        except OSError as e:
                if e.errno == errno.EEXIST:
                    logging.error("\nERROR: Folders already exist, won't "
                                  "overwrite.")
                    sys.exit()
                else:
                    raise

    if cart_coords:
        logging.info('Writing band structure in Cartesian coordinates...\n'
                     'Remember to run full-potential calc with --band and '
                     'NOT --band~mq')
    else:
        logging.info('Writing band structure in direct coordinates...\n'
                     'Remember to run full-potential calc with --band~mq.')
    if labels is None:
        with open(os.path.join(path, 'syml.' + ext), 'w') as f:
            for kpt in kpoints:
                f.write('{0:11.8f} {1:11.8f} {2:11.8f}\n'.format(kpt[0],
                                                                 kpt[1],
                                                                 kpt[2]))
    else:
        label_positions = [i for i, l in enumerate(labels) if l != '']
        special_points = [kpoints[i] for i in label_positions]
        segment_samples = [label_positions[i + 1] - label_positions[i] + 1
                           for i in range(len(label_positions) - 1)]
        with open(os.path.join(path, 'syml.' + ext), 'w') as f:
            for i, samples in enumerate(segment_samples):
                if samples == 2:
                    continue   # Don't add segments between branches
                f.write('{0:5d}    {1:11.8f} {2:11.8f} {3:11.8f}    '
                        '{4:11.8f} {5:11.8f} {6:11.8f}    {7} to {8}\n'.format(
                            samples,
                            special_points[i][0], special_points[i][1],
                            special_points[i][2],
                            special_points[i + 1][0], special_points[i + 1][1],
                            special_points[i + 1][2],
                            labels[label_positions[i]],
                            labels[label_positions[i + 1]]))
            f.write('    0 0 0 0 0 0 0\n')


def labels_from_syml(syml_file):
    """Read special point locations from Questaal syml file

    Each line of file should have format:

        NPTS  X1 Y1 Z1  X2 Y2 Z2  label1 to label2

    except the last line which is seven 0 characters. Whitespace amounts
    are flexible. If a label is defined more than once, the last definition
    will be used.

    Args:
        syml_file (:obj:`str`): Path to questaal syml.ext input file. The
            locations and labels of special points are read from this file.

    """
    labels = {}

    with open(syml_file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        npts, x1, y1, z1, x2, y2, z2, *label_text = line.split()
        if len(label_text) < 3:
            pass
        else:
            kpt1 = tuple(map(float, (x1, y1, z1)))
            kpt2 = tuple(map(float, (x2, y2, z2)))

            label_text = ' '.join(label_text)  # Undo previous split
            label1_label2 = label_text.split(' to ')
            if len(label1_label2) != 2:
                raise ValueError("Not clear how to interpret labels from "
                                 "this line: {}".format(line))
            label1, label2 = label1_label2
            labels.update({label1: kpt1, label2: kpt2})
    return labels


def read_dos(pdos_file=None, tdos_file=None, site_file=None,
             ry=True, total_only=False, gaussian=None,
             elements=None, lm_orbitals=None, atoms=None):
    """Read DOS file data

    Construct Pymatgen Dos objects for total dos and orbital contributions as
    available from Questaal dos.ext file

    Args:
        pdos_file (:obj:`str`):
            Path to DOS file from Questaal, as generated by LMF with ``--dos``
            option (to obtain total DOS) or with LMDOS from moms.ext (to obtain
            PDOS). It is assumed in either case that the default mode=2 was
            used (i.e. DOS is resolved by site, l and m). If there is only one
            channel of data, this is interpreted as a total DOS.
        tdos_file (:obj:`str`):
            Path to TDOS file from Questaal, as generated by LMF with ``--dos``
            option BEFORE LMDOS is used to generate a PDOS from *moms.ext*.
        site_file (:obj:`str`, optional):
            Path to site.ext file from Questaal. This is used to identify the
            sites and obtain structure information for PDOS manipulation.
        ry (:obj:`bool`, optional):
            Convert energy values from Ry to eV
        elements (:obj:`dict`, optional): The elements and orbitals to extract
            from the projected density of states. Should be provided as a
            :obj:`dict` with the keys as the element names and corresponding
            values as a :obj:`tuple` of orbitals. For example, the following
            would extract the Bi s, px, py and d orbitals::

                {'Bi': ('s', 'px', 'py', 'd')}

            If an element is included with an empty :obj:`tuple`, all orbitals
            for that species will be extracted. If ``elements`` is not set or
            set to ``None``, all elements for all species will be extracted.
        lm_orbitals (:obj:`dict`, optional): The orbitals to decompose into
            their lm contributions (e.g. p -> px, py, pz). Should be provided
            as a :obj:`dict`, with the elements names as keys and a
            :obj:`tuple` of orbitals as the corresponding values. For example,
            the following would be used to decompose the oxygen p and d
            orbitals::

                {'O': ('p', 'd')}

        atoms (:obj:`dict`, optional): Which atomic sites to use when
            calculating the projected density of states. Should be provided as
            a :obj:`dict`, with the element names as keys and a :obj:`tuple` of
            :obj:`int` specifying the atomic indices as the corresponding
            values. The elemental projected density of states will be summed
            only over the atom indices specified. If an element is included
            with an empty :obj:`tuple`, then all sites for that element will
            be included. The indices are 0 based for each element specified in
            the site.ext. For example, the following will calculate the density
            of states for the first 4 Sn atoms and all O atoms in the
            structure::

                {'Sn': (1, 2, 3, 4), 'O': (, )}

            If ``atoms`` is not set or set to ``None`` then all atomic sites
            for all elements will be considered.
        gaussian (:obj:`float`, optional): Broaden the density of states using
            convolution with a gaussian function. This parameter controls the
            sigma or standard deviation of the gaussian distribution.
        total_only (:obj:`bool`, optional): Only extract the total density of
            states. Defaults to ``False``.

    Returns:
        (:obj:`tuple`):
            *dos*, *pdos*: *dos* is a
            :obj:`~pymatgen.electronic_structure.dos.Dos` object containing the
            total density of states and ``pdos`` is a :obj:`dict` of
            :obj:`dict` mapping the elements and their orbitals to
            :obj:`~pymatgen.electronic_structure.dos.Dos` objects. For
            example::

                 {
                     'Bi': {'s': Dos, 'p': Dos ... },
                     'S': {'s': Dos}
                 }
    """

    def _read_states(f, nlines):
        lines = [f.readline() for i in range(nlines)]
        # This statement is very "functional programming"; read it
        # backwards.  List of split lines is "flattened" by chain into
        # iterator of values; this is fed into map to make floats and
        # stored to a list
        return list(map(float, chain(*(line.split() for line in lines))))

    def _read_dos_data(filename, check_tdos=False):
        with zopen(filename, 'rt') as f:
            (emin, emax, ne, nchan,
             nsp, efermi, delta, fmt) = f.readline().split()

            if check_tdos and int(nchan) != 1:
                raise ValueError('File {} is not a TDOS: contains {} channels'
                                 ''.format(filename, nchan))

            ne, nchan, nsp, fmt = map(int, (ne, nchan, nsp, fmt))
            emin, emax, efermi, delta = map(float, (emin, emax, efermi, delta))
            if ry:
                emin, emax, efermi, delta = [x * _ry_to_ev for x in
                                             (emin, emax, efermi, delta)]

            nlines = ceil(ne / 5)
            data = [_read_states(f, nlines) for _ in range(nsp * nchan)]

        energies = np.linspace(emin, emax, ne)

        return energies, data, efermi, nsp

    def _get_tdos(filename):
        energies, data, efermi, nsp = _read_dos_data(tdos_file,
                                                     check_tdos=True)

        if nsp == 1:
            # densities = np.array(data[0]).reshape((ne, 1))
            densities = {Spin.up: data[0]}
        elif nsp == 2:
            #densities = np.concatenate((data[0], data[1]), axis=1)
            densities = {Spin.up: data[0], Spin.down: data[1]}
        else:
            raise ValueError('There can\'t be {} spin channels, that makes '
                             'no sense!'.format(nsp))
        return Dos(efermi, energies, densities)

    def _get_cdos(pdos_file, site_file):
        if not (pdos_file and site_file):
            raise ValueError('Both dos.ext and site.ext are needed for PDOS')
        energies, data, efermi, nsp = _read_dos_data(pdos_file)
        if nsp == 1:
            fake_tdos = Dos(efermi, energies, {Spin.up: 0 * energies})
            spins = (Spin.up,)
        elif nsp == 2:
            fake_tdos = Dos(efermi, energies, {Spin.up: 0 * energies,
                                               Spin.down: 0 * energies})
            spins = (Spin.up, Spin.down)
        else:
            raise ValueError('There can\'t be {} spin channels, that makes '
                             'no sense!'.format(nsp))

        site_data = QuestaalSite.from_file(site_file)
        structure = site_data.structure

        pdoss = {}
        for site, orbital, spin in product(range(len(site_data.sites)),
                                           range(16),  # forget about g orbs
                                           range(len(spins))):
            if structure.sites[site] not in pdoss:
                pdoss.update({structure.sites[site]: {}})
            if Orbital(orbital) not in pdoss[structure.sites[site]]:
                pdoss[structure.sites[site]].update({Orbital(orbital): {}})

            pdoss[structure.sites[site]][Orbital(orbital)].update(
                {spins[spin]: data[site * (25 * nsp) + orbital * nsp + spin]})

        return CompleteDos(structure, fake_tdos, pdoss)

    if not (tdos_file or pdos_file):
        raise ValueError('Need some DOS data: provide a dos.ext')

    elif tdos_file and not pdos_file:
        tdos = _get_tdos(tdos_file)
        cdos = None

    elif pdos_file and not tdos_file:
        with zopen(pdos_file, 'rt') as f:
            nchan = int(f.readline().split()[3])
        if nchan == 1:
            logging.info("No explicit TDOS and file {} contains 1 channel: "
                         "this is a TDOS.".format(pdos_file))

            tdos = _get_tdos(pdos_file)
            cdos = None
        else:
            tdos = None
            cdos = _get_cdos(pdos_file, site_file)

    else:
        tdos = _get_tdos(tdos_file)
        cdos = _get_cdos(pdos_file, site_file)

    if gaussian:
        if tdos is not None:
            tdos.densities = tdos.get_smeared_densities(gaussian)
        if cdos is not None:
            for site in cdos.pdos:
                for orbital in cdos.pdos[site]:
                    cdos.pdos[site][orbital] = cdos.get_site_orbital_dos(
                        site, orbital).get_smeared_densities(gaussian)

    # Combine data into selected orbitals: default is by species type
    if cdos is None:
        pdos = {}
    else:
        pdos = sumo.electronic_structure.dos.get_pdos(
            cdos, lm_orbitals=lm_orbitals, atoms=atoms, elements=elements)

    return (tdos, pdos)


def band_structure(bnds_file, lattice, labels=None, alat=1,
                   coords_are_cartesian=False):
    """Read band structure data

    Args:
        bands_file (:obj:`str`): Path to questaal bnds.ext output file. The
            k-point positions and eigenvalues are read from this file.
        labels (:obj:`dict`): Dict of special point locations and labels. This
            is generally obtained from a syml.ext file using
            :obj:`sumo.io.questaal.labels_from_syml`.
        alat (:obj:`float`): Lattice scaling parameter defined in site.ext
        lattice (:obj:`pymatgen.core.lattice.Lattice`)
        coords_are_cartesian (:obj:`bool`): bnds.ext file is in Cartesian
            coordinates.

    Returns:
        :obj:`pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`
    """

    # Implementation note: pymatgen band structure expects an eigenvalue array
    # ordered [band][kpt] but the data will be read in one kpt at a time as
    # this is the format of the bnd data file.
    #
    # We will first build a nested list [[bnd1_kpt1, bnd2_kpt1, ...],
    #                                    [bnd1_kpt2, bnd2_kpt2, ...], ...]
    # then convert to a numpy array (i.e. with each row for a kpoint)
    # and then transpose the array to obtain desired formats

    with zopen(bnds_file, 'r') as f:
        kpoints = []

        # Read heading, get metadata and check no orbital projections used
        nbands, efermi, n_color_wts, *_ = f.readline().split()

        if int(n_color_wts) > 0:
            raise NotImplementedError("Band data includes orbital data: "
                                      "this format is not currently supported."
                                      )

        nbands, efermi = int(nbands), float(efermi)
        eig_lines = ceil(nbands / 10)

        # Check if first two kpts are the same: if so, assume there are two
        # spin channels
        _ = f.readline()
        kpt1 = list(map(float, f.readline().split()))

        for line in range(eig_lines):      # Skip over the eigenvalues
            _ = f.readline()                     # for now: re-read file later
        kpt2 = list(map(float, f.readline().split()))
        if len(kpt1) != 3 or len(kpt2) != 3:
            raise AssertionError()

        if kpt1 == kpt2:
            spin_pol = True
        else:
            spin_pol = False

    logging.info("Reading Questaal band structure header...")
    logging.info("nbands: {}, efermi / Ry: {}, spin channels: {}".format(
        nbands, efermi, (2 if spin_pol else 1)))
    logging.info("Reading band structure eigenvalues...")

    # Key info has been collected: re-read file for kpts and eigenvalues
    def _read_eigenvals(f, nlines):
        lines = [f.readline() for i in range(nlines)]
        # This statement is very "functional programming"; read it
        # backwards.  List of split lines is "flattened" by chain into
        # iterator of values; this is fed into map to make floats and
        # stored to a list
        return list(map(float, chain(*(line.split() for line in lines))))

    with zopen(bnds_file, 'r') as f:
        _ = f.readline()
        if spin_pol:  # Need to read two spin channels
            block_nkpts = int(f.readline().strip()) // 2
            eigenvals = {Spin.up: [], Spin.down: []}
        else:
            block_nkpts = int(f.readline().strip())
            eigenvals = {Spin.up: []}

        while block_nkpts > 0:  # File should be terminated with a 0
            for i in range(block_nkpts):
                kpoint = list(map(float, f.readline().split()))
                kpoints.append(np.array(kpoint) / (alat * _bohr_to_angstrom))

                eigenvals[Spin.up].append(_read_eigenvals(f, eig_lines))

                if spin_pol:
                    spin_down_kpoint = list(map(float, f.readline().split()))
                    if spin_down_kpoint != kpoint:
                        raise AssertionError(
                            "File interpreted as spin-polarised, but this"
                            " kpoint only has one entry: {}".format(kpoint))
                    eigenvals[Spin.down].append(_read_eigenvals(f, eig_lines))

            block_nkpts = int(f.readline().strip())
            if spin_pol:
                block_nkpts = block_nkpts // 2

    # Transpose matrix to arrange by band and convert to eV from Ry
    eigenvals = {key: np.array(data).T * _ry_to_ev
                 for key, data in eigenvals.items()}
    efermi *= _ry_to_ev

    # Convert labels to Cartesian coordinates
    # [a, b, c] [ax, ay, az] = [a.ax + b.bx + x.cx, a.ay + b.by,...] = [x y z]
    #           |bx, by, bz|
    #           [cx, cy, cz]

    labels = labels or {} # Initialise a dict if null value
                          # (Avoids oddness with mutable option defaults)
    if coords_are_cartesian:
        logging.info("Cartesian coordinates, scaling by ALAT = {}".format(alat))
        for label, coords in labels.items():
            labels[label] = np.array(coords) / (alat * _bohr_to_angstrom)
    else:
        for label, coords in labels.items():
            labels[label] = np.dot(
                coords, lattice.reciprocal_lattice_crystallographic.matrix)

    # Data in bnds file seems to always be Cartesian
    return BandStructureSymmLine(kpoints, eigenvals,
                                 lattice.reciprocal_lattice_crystallographic,
                                 efermi, labels,
                                 coords_are_cartesian=True)


def dielectric_from_file(filename, out_filename=None):
    """Detect Questaal optics file type and dispatch to appropriate reader

    Args:
        filename (:obj:`str`):
            Path to *opt.ext* or *eps_BSE.out* data file from Questaal. If the
            filename contains "eps_BSE" the *eps_BSE.out* format is assumed,
            otherwise the file is treated as an *opt.ext*.
        out_filename (:obj:`str`, optional):
            Filename for writing out the calculated dielectric data. This is
            only used for *opt.ext* inputs; energy values are converted to eV
            and a Real component is obtained by the Kramers-Kronig relation.
    """

    if 'eps_BSE' in filename:
        return dielectric_from_BSE(filename)
    else:
        return dielectric_from_opt(filename, out_filename=out_filename)


def dielectric_from_BSE(filename):
    """Read a Questaal eps_BSE.out file and return dielectric function

    eps_BSE files only provide a scalar complex number; this is converted to a
    diagonal matrix for compatibility with other routines. Unlike the opt.exe
    however the real part is provided so we do not need to make a
    Kramers-Kronig transformation.

    Args:
        filename (:obj:`float`):
            Path to ``eps_BSE.out`` output of bethesalpeter calculation.

    Returns:
        :obj:`tuple`
            The format imitates the ``dielectric`` attribute of
            :obj:`pymatgen.io.vasp.outputs.Vasprun`: a tuple of the form::

                ([energy1, energy2, ...],
                 [[real_xx_1, real_yy_1, real_zz_1,
                   real_xy_1, real_yz_1, real_xz_1],
                  [real_xx_2, real_yy_2, real_zz_2,
                   real_xy_2, real_yz_2, real_xz_2], ...],
                 [[imag_xx_1, imag_yy_1, imag_zz_1,
                   imag_xy_1, imag_yz_1, imag_xz_1],
                  [imag_xx_2, imag_yy_2, imag_zz_2,
                   imag_xy_2, imag_yz_2, imag_xz_2], ...])

            As only a scalar is given in eps_BSE, the vectors will in fact
            be e.g.::

            [real_1, real_1, real_1, 0, 0, 0]
    """
    data = np.genfromtxt(filename, comments='#')
    if data.shape[1] == 3:
        pass
    else:
        raise ValueError("Not sure how to interpret {}; expected "
                         "3 columns. "
                         "If this isn't an eps_BSE file, please don't put"
                         "eps_BSE in the filename!".format(filename))

    energy = list(data[:, 0].T)
    real = [[r, r, r, 0, 0, 0] for r in data[:, 1]]
    imag = [[i, i, i, 0, 0, 0] for i in data[:, 2]]

    return(energy, real, imag)


def dielectric_from_opt(filename, cshift=1e-6, out_filename=None):
    """Read a Questaal opt.ext file and return dielectric function

    opt.ext files only provide x, y, z components so the off-diagonal terms are
    set to zero.

    Args:
        filename (:obj:`float`):
            Path to ``opt.ext`` output of LMF optics calculation. If this is
            split into spin channels, the values will be recombined to a single
            (observable) channel.
        cshift (:obj:`float`, optional):
            A small imaginary element used in Kramers-Kronig integration.
        out_filename (:obj:`float`, optional):
            Path to write tabulated dielectric data with columns::

                energy(eV) real_xx real_yy real_zz imag_xx imag_yy imag_zz

            If None, no file is written.

    Returns:
        :obj:`tuple`
            The format imitates the ``dielectric`` attribute of
            :obj:`pymatgen.io.vasp.outputs.Vasprun`: a tuple of the form::

                ([energy1, energy2, ...],
                 [[real_xx_1, real_yy_1, real_zz_1,
                   real_xy_1, real_yz_1, real_xz_1],
                  [real_xx_2, real_yy_2, real_zz_2,
                   real_xy_2, real_yz_2, real_xz_2], ...],
                 [[imag_xx_1, imag_yy_1, imag_zz_1,
                   imag_xy_1, imag_yz_1, imag_xz_1],
                  [imag_xx_2, imag_yy_2, imag_zz_2,
                   imag_xy_2, imag_yz_2, imag_xz_2], ...])

            The off-diagonal (xy, yz, xz) terms are not given in ``opt.ext``
            and are set to zero. Only the imaginary part is provided so the
            real part is constructed by a Kramers-Kronig trasformation.
            Energy units are converted from Ry to eV.
    """

    data = np.genfromtxt(filename, skip_header=1)
    if data.shape[1] == 4:
        pass
    elif data.shape[1] == 7:
        data = np.hstack(data[:, :1], data[1:4] + data[4:7])
    else:
        raise ValueError("Not sure how to interpret {}; expected "
                         "4 or 7 columns.".format(filename))

    data[:, 0] *= _ry_to_ev
    de = data[1, 0] - data[0, 0]
    eps_imag = [[[row[1], 0, 0], [0, row[2], 0], [0, 0, row[3]]]
                for row in data]
    eps_real = kkr(de, eps_imag)

    # Re-shape to XX YY ZZ XY YZ XZ format
    eps_imag = np.array(eps_imag).reshape(len(eps_imag), 9)
    eps_imag = eps_imag[:, [0, 4, 8, 1, 5, 2]]
    eps_real = eps_real.reshape(len(eps_real), 9)
    eps_real = eps_real[:, [0, 4, 8, 1, 5, 2]]

    if out_filename is not None:
        with open(out_filename, 'wt') as f:
            for e, r, i in zip(data[:, 0].flatten(), eps_real, eps_imag):
                f.write((' '.join(['{:10.8f}'] * 7)).format(e, *r[:3], *i[:3]))
                f.write('\n')

    return (data[:, 0].flatten(), eps_real, eps_imag)
