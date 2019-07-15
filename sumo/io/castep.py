import collections
import logging
import os
import re
from monty.io import zopen
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

_bohr_to_angstrom = 0.5291772
_ry_to_ev = 13.605693009

Tag = collections.namedtuple('Tag', 'value comment')
Block = collections.namedtuple('Block', 'values comments')

class CastepCell(object):
    """Structure information and more: CASTEP seedname.cell file

    Usually this will be instantiated with the
    :obj:`sumo.io.castep.CastepCell.from_file()` method

    Tags are stored as a dict, where each key is a tag and each value is a
    Tag NamedTuple with the attributes 'value' and 'comment'.

    Blocks are also stored as a dict, where each key is the block name
    (e.g. 'bs_kpoint_list') and the value is a Block NamedTuple with attributes
    'values' (a list of lists) and 'comments' (a commensurate list of str).

    """

    def __init__(self, blocks={}, tags={}):
        self.blocks = blocks
        self.tags = tags

    @classmethod
    def from_file(cls, filename):

        with zopen(filename, 'rt') as f:
            lines = [line.strip() for line in f]

        # Remove lines which are entirely commented or empty
        def _is_not_empty_line(line):
            if len(line) == 0:
                return False
            elif line[0] in '#;!':
                return False
            elif len(line) > 6 and line[:7] == 'COMMENT':
                return False
            else:
                return True

        lines = list(filter(_is_not_empty_line, lines))

        tags, blocks = {}, {}
        current_block_values, current_block_comments = [], []
        in_block = False
        current_block_label = ''

        for line in lines:
            if len(line.split()) == 0:
                continue
            elif line[:6].lower() == '%block':
                if in_block:
                    raise IOError('Cell file contains nested blocks. '
                                  'This possibility was not anticipated.')
                else:
                    current_block_label = line.split()[1].lower()
                    in_block = True
                    continue
            elif line[:9].lower() == '%endblock':
                if line.split()[1].lower() != current_block_label:
                    raise IOError('Endblock {} does not match current block '
                                  '{}: cannot interpret cell file.'.format(
                                      line.split()[1], current_block_label))
                if in_block:
                    blocks[current_block_label] = Block(current_block_values,
                                                        current_block_comments)

                    current_block_values, current_block_comments = [], []
                    in_block = False
                else:
                    raise IOError('Cannot cope with line {}: not currently in '
                                  'a block.'.format(line))
            elif in_block:
                comment_split_line = re.split('[#;!]', line)
                if len(comment_split_line) == 1:
                    current_block_values.append(line.split())
                    current_block_comments.append('')
                else:
                    values = comment_split_line[0]
                    current_block_values.append(values.split())
                    current_block_comments.append(line[len(values):])
            else:
                comment_split_line = re.split('[#;!]', line)
                if len(comment_split_line) == 1:
                    comment = ''
                    line = line.split()
                else:
                    comment = line[len(comment_split_line[0]):]
                    line = comment_split_line[0].split()
                    line = [item for item in line if item not in ':=']

                if len(line) == 1:
                    data = ''
                else: data = line[1:]

                tag = line[0].lower()
                if tag[-1] in ':=':
                    tag = tag[:-1]

                tags[tag] = Tag(data, comment)

        return cls(tags=tags, blocks=blocks)


def band_structure(bands_file, cell_file=None):
    """Convert band structure data from CASTEP to Pymatgen/Sumo format

    Args:
        bands_file (:obj:`str`): Path to CASTEP prefix.bands output file. The
            lattice parameters, k-point positions and eigenvalues are read from
            this file.
        cell_file (:obj:`str`, optional): Path to CASTEP prefix.cell input
            file. If found, the positions of high-symmetry points are read from
            the ``bs_kpoint_path`` block in this file.

    Returns:
        :obj:`pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`

    """

    logging.info("Reading CASTEP band structure header...")
    header = read_bands_header(bands_file)

    if header['n_spins'] == 1:
        logging.info("nbands: {}, efermi / Ry: {}, spin channels: 1".format(
            header['n_bands'][0],
            header['e_fermi'][0]))
    elif header['n_spins'] == 2:
        logging.info("nbands: {}, efermi / Ry: {}, spin channels: 2".format(
            header['n_bands'],
            header['e_fermi']))
    else:
        raise ValueError('Should only be 1 or 2 spin channels!')

    if len(header['e_fermi']) == 2:
        if abs(header['e_fermi'][1] - header['e_fermi'][1]) > 1e-8:
            raise NotImplementedError("Different Fermi energy in each channel."
                                      " I have no idea how to handle that.")

    lattice = Lattice(header['lattice_vectors'])
    lattice = Lattice(lattice.matrix * _bohr_to_angstrom)

    logging.info("Reading band structure eigenvalues...")
    kpoints, eigenvalues = read_bands_eigenvalues(bands_file, header)

    if cell_file is  None:
        labels = {}
    else:
        labels = labels_from_cell(cell_file)

    return BandStructureSymmLine(kpoints, eigenvalues,
                                 lattice.reciprocal_lattice_crystallographic,
                                 header['e_fermi'][0] * _ry_to_ev * 2,
                                 labels,
                                 coords_are_cartesian=False)


def labels_from_cell(cell_file):
    """Read special k-point positions/labels from CASTEP cell file

    Args:
        cell_file (:obj:`str`): Path to CASTEP prefix.cell input file
            high-symmetry points will be read from a block beginning with

                %block bs_kpoint_path

            and ending with

                %endblock bs_kpoint path

            where each special-point line is formatted

                k1 k2 k3  ! label

            where k1, k2, k3 are given in fractional coordinates of the
            reciprocal lattice.

    Returns:
        {:obj:`str`: :obj:`tuple`, ...}

    """

    labels = {}

    blockstart = re.compile('^%block\s+bs_kpoint_(path|list)')
    blockend = re.compile('^%endblock\s+bs_kpoint_(path|list)')

    with zopen(cell_file, 'r') as f:
        line = ''
        while blockstart.match(line.lower()) is None:
            line = f.readline()

        line = f.readline()  # Skip past block start line
        while blockend.match(line.lower()) is None:
            kpt = tuple(map(float, line.split()[:3]))
            if len(line.split()) > 3:
                label = line.split()[-1]
                if label.lower() in ('g', 'gamma'):
                    label = '\Gamma'
                labels[label] = kpt
            line = f.readline()

    return labels

def read_bands_header(bands_file):
    """Read CASTEP bands file header, get lattice vectors and metadata

    Args:
        bands_file (:obj:`str`): Path to CASTEP prefix.bands output file

    Returns:
        :obj:`dict`
            CASTEP bands metadata of the form::

            {'lattice': [[ax, ay, az], [bx, by, bz], [cx, cy, cz]],
             'n_kpoints': n_kpoints,
             'n_spins': n_spins,
             'n_electrons': [n_spin1, n_spin2],
             'n_bands': [n_bands_spin1, n_bands_spin2],
             'e_fermi': [e_fermi_spin1, e_fermi_spin2]}

    """

    header = {}
    with zopen(bands_file, 'r') as f:
        # Read in header information a line at a time

        header['n_kpoints'] = int(f.readline().split()[-1])
        header['n_spins'] = int(f.readline().split()[-1])
        header['n_electrons'] = [float(x) for x in f.readline().split()[3:]]
        header['n_bands'] = [int(x) for x in f.readline().split()[3:]]

        fermi_line = f.readline().split()
        if fermi_line[3] != 'atomic':
            raise NotImplementedError('Band data not in atomic units')
        header['e_fermi'] = [float(x) for x in fermi_line[5:]]

        vectors_header = f.readline().strip()
        if vectors_header != 'Unit cell vectors':
            raise AssertionError('CASTEP bands file not formatted as expected', vectors_header)

        lattice_vectors = [[float(x) for x in f.readline().split()]
                           for _ in range(3)]
        for row in lattice_vectors:
            if len(row) != 3:
                raise AssertionError('Unit cell vectors not read correctly')
        header['lattice_vectors'] = lattice_vectors

    return header


def read_bands_eigenvalues(bands_file, header):
    """Read eigenvalue data from a CASTEP bands file

    Args:
        bands_file (:obj:`str`): Path to CASTEP prefix.bands output file
        header (:obj:`dict`): Band metadata, as obtained by reading the bands
            file with :obj:`read_bands_header`.

    Returns:
        (:obj:`list`, :obj:`np.ndarray`)
            kpoints list and eigenvalue array ordered by [band][kpt]. Each
            k-point is expressed as numpy array of length 3. Eigenvalues are
            converted from Ry to eV.

    """

    # Implementation note: pymatgen band structure expects an eigenvalue array
    # ordered [band][kpt] but the data will be read in one kpt at a time as
    # this is the format of the bnd data file.
    #
    # We will build a nested list [[bnd1_kpt1, bnd2_kpt1, ...],
    #                              [bnd1_kpt2, bnd2_kpt2, ...], ...]
    # then convert to a numpy array (i.e. with each row for a kpoint)
    # and then transpose the array to obtain desired formats
    #
    # CASTEP can write the k-points out-of-order, so we need to keep track of
    # the k-point indices and re-order the lists before array conversion.

    kpoints = []
    if header['n_spins'] == 1:
        eigenvals = {Spin.up: []}
    else:
        eigenvals = {Spin.up: [], Spin.down: []}

    with zopen(bands_file, 'r') as f:
        for i in range(9):
            _ = f.readline()  # Skip header

        for i_kpt in range(header['n_kpoints']):
            # The first "k-point" column is actually the sorting index
            kpoints.append(np.array(
                [float(x) for x in f.readline().split()[1:5]]))
            _ = f.readline() # Skip past "Spin component 1"
            spin1_bands = ([kpoints[-1][0]] +     # Sorting key goes in col 0
                           [float(f.readline())
                            for _ in range(header['n_bands'][0])])
            eigenvals[Spin.up].append(spin1_bands)

            if header['n_spins'] == 2:
                _ = f.readline() # Skip past "Spin component 2"
                spin2_bands = ([kpoints[-1][0]] + # Sorting key goes in col 0
                               [float(f.readline())
                                for _ in range(header['n_bands'][0])])
                eigenvals[Spin.down].append(spin2_bands)

    # Sort kpoints and trim off sort keys
    kpoints = sorted(kpoints, key=(lambda k: k[0]))
    kpoints = [k[1:] for k in kpoints]

    # Sort eigenvalues and trim off sort keys
    for key, data in eigenvals.items():
        data = np.array(data)
        eigenvals[key] = data[data[:, 0].argsort(), 1:]
    # Transpose matrix to arrange by band and convert to eV from Ha
    eigenvals = {key: data.T * _ry_to_ev * 2
                 for key, data in eigenvals.items()}

    return kpoints, eigenvals

def write_kpoint_files(filename, kpoints, labels, make_folders=False,
                       kpts_per_split=None, directory=None):
    r"""Write the k-points data to files.

    Folders are named as 'split-01', 'split-02', etc ...
    .cell files are named band_split_01.cell etc ...

    Cartesian coordinates are not supported by CASTEP for band paths.

    Unlike VASP, CASTEP band structure calculations can deal with two
    separate k-point meshes: a mesh for the SCF convergence and a series of
    points (or path segments) for the band structure calculation. It is assumed
    that the existing .cell file contains a specification of the SCF k-point
    mesh, and Sumo will append the band path information.

    Args:
        filename (:obj:`str`): Path to CASTEP structure (.cell) file.
        kpoints (:obj:`numpy.ndarray`): The k-point coordinates along the
            high-symmetry path. For example::

                [[0, 0, 0], [0.25, 0, 0], [0.5, 0, 0], [0.5, 0, 0.25],
                [0.5, 0, 0.5]]

        labels (:obj:`list`) The high symmetry labels for each k-point (will be
            an empty :obj:`str` if the k-point has no label). For example::

                ['\Gamma', '', 'X', '', 'Y']

        make_folders (:obj:`bool`, optional): Generate folders and copy in
            param file from the current directory, setting task=BandStructure.

        kpts_per_split (:obj:`int`, optional): If set, the k-points are split
            into separate k-point files (or folders) each containing the number
            of k-points specified. This is useful for hybrid band structure
            calculations where it is often intractable to calculate all
            k-points in the same calculation.

        directory (:obj:`str`, optional): The output file directory.
        """

    if kpts_per_split:
        kpt_splits = [kpoints[i:i+kpts_per_split] for
                      i in range(0, len(kpoints), kpts_per_split)]
        label_splits = [labels[i:i+kpts_per_split] for
                        i in range(0, len(labels), kpts_per_split)]
    else:
        kpt_splits = [kpoints]
        label_splits = [labels]

    kpt_txt_blocks = []
    for kpt_split, label_split in zip(kpt_splits, label_splits):
        # kpt weights don't actually matter for Castep band structure;
        # Create an array of zeroes so Pymatgen has something to work with
        kpt_weights = np.zeros(len(kpt_split))

        txt_block = '\n'.join(
            ['%block bs_kpoint_list']
            + ['{10.8f} {10.8f} {10.8f} ! {}'.format(k[0], k[1], k[2], l)
               for k, l in zip(kpt_split, label_split)]
            + '%endblock bs_kpoint_list\n\n')

        kpt_txt_blocks.append(kpt_txt_blocks)

    pad = int(math.floor(math.log10(len(kpt_files)))) + 2
    if make_folders:
        # Derive param file name by changing extension from cell file
        param_file = '.'.join(filename.split('.')[:-1] + ['param'])

        with open(filename, 'r') as f:
            cell_file_txt = f.readlines()

        for i, kpt_txt_block in enumerate(kpt_txt_blocks):
            split_folder = 'split-{}'.format(str(i + 1).zfill(pad))
            if directory:
                split_folder = os.path.join(directory, split_folder)

            copy_param(param_file, split_folder, task='BandStructure')

            with open(os.path.join(split_folder,
                                   os.path.basename(filename)), 'w') as f:
                f.write(cell_file_txt)
                f.write('\n')
                f.write(kpt_txt_block)

    else:
        for i, kpt_txt_block in enumerate(kpt_txt_blocks):
            if len(kpt_txt_blocks) > 1:
                cell_filename = 'band_split_{:0d}.cell'.format(i + 1)
            else:
                cell_filename = 'band.cell'
            if directory:
                cell_filename = os.path.join(directory, cell_filename)
            with open(cell_filename, 'w') as f:
                f.write(kpt_txt_block)

def copy_param(filename, folder, task=None):
    """Copy CASTEP .param file, modifying Task if needed

    Args:
        filename (:obj:`str`): Path to CASTEP .param file

    """

    output_filename = os.path.join(folder, os.path.basename)
    if os.path.isfile(filename) and task is None:
        shutil.copyfile(filename, os.path.join(folder, filename))

    elif os.path.isfile(filename):
        with open(filename, 'r')  as infile:
            with open(os.path.join(folder, output_filename), 'w') as outfile:
                for line in infile:
                    if line.strip().lower()[:4] == 'task':
                        outfile.write('Task : BandStructure\n')
                    else:
                        outfile.write(line)
    else:
        logging.warning("Cannot find param file {}, skipping".format(filename))
