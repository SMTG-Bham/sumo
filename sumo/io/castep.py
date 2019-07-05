import logging
import re
from monty.io import zopen
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine

_bohr_to_angstrom = 0.5291772
_ry_to_ev = 13.605693009

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
