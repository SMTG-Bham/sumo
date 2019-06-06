# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
A script to plot phonon band structure diagrams.

TODO:
 * automatically plot dos if present in band.yaml
 * primitive axis determination (try symmetrise->spglib->PA then apply
   transform on original cell and see if it works.
 * make band structure from vasprun displacement/DFPT files
 * deal with magnetic moments
 * Read FORCE_CONSANTS or force_constants.hdf5
 * change frequency unit
 * read settings from phonopy config file
 * prefix file names
"""

import os
from os.path import isfile
import sys
import logging
import argparse
import numpy as np

import warnings
warnings.filterwarnings("ignore", category=FutureWarning,
                        module="h5py")

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rcParams

from phonopy.units import VaspToTHz, VaspToEv, VaspToCm

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.phonopy import get_ph_bs_symm_line

from sumo.phonon.phonopy import load_phonopy
from sumo.symmetry.kpoints import get_path_data
from sumo.plotting.phonon_bs_plotter import SPhononBSPlotter


__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Jan 17, 2018"


def phonon_bandplot(filename, poscar=None, prefix=None, directory=None,
                    dim=None, born=None, qmesh=None, spg=None,
                    primitive_axis=None, line_density=60, units='THz',
                    symprec=0.01, mode='bradcrack', kpt_list=None,
                    eigenvectors=False, labels=None, height=6., width=6.,
                    style=None, no_base_style=False,
                    ymin=None, ymax=None, image_format='pdf', dpi=400,
                    plt=None, fonts=None, dos=None):
    """A script to plot phonon band structure diagrams.

    Args:
        filename (str): Path to phonopy output. Can be a band structure yaml
            file, ``FORCE_SETS``, ``FORCE_CONSTANTS``, or
            ``force_constants.hdf5``.
        poscar (:obj:`str`, optional): Path to POSCAR file of unitcell. Not
            required if plotting the phonon band structure from a yaml file. If
            not specified, the script will search for a POSCAR file in the
            current directory.
        prefix (:obj:`str`, optional): Prefix for file names.
        directory (:obj:`str`, optional): The directory in which to save files.
        born (:obj:`str`, optional): Path to file containing Born effective
            charges. Should be in the same format as the file produced by the
            ``phonopy-vasp-born`` script provided by phonopy.
        qmesh (:obj:`list` of :obj:`int`, optional): Q-point mesh to use for
            calculating the density of state. Formatted as a 3x1 :obj:`list` of
            :obj:`int`.
        spg (:obj:`str` or :obj:`int`, optional): The space group international
            number or symbol to override the symmetry determined by spglib.
            This is not recommended and only provided for testing purposes.
            This option will only take effect when ``mode = 'bradcrack'``.
        primitive_matrix (:obj:`list`, optional): The transformation matrix
            from the conventional to primitive cell. Only required when the
            conventional cell was used as the starting structure. Should be
            provided as a 3x3 :obj:`list` of :obj:`float`.
        line_density (:obj:`int`, optional): Density of k-points along the
            path.
        units (:obj:`str`, optional): Units of phonon frequency. Accepted
            (case-insensitive) values are Thz, cm-1, eV, meV.
        symprec (:obj:`float`, optional): Tolerance for space-group-finding
            operations
        mode (:obj:`str`, optional): Method used for calculating the
            high-symmetry path. The options are:

            bradcrack
                Use the paths from Bradley and Cracknell. See [brad]_.

            pymatgen
                Use the paths from pymatgen. See [curt]_.

            seekpath
                Use the paths from SeeK-path. See [seek]_.

        kpt_list (:obj:`list`, optional): List of k-points to use, formatted as
            a list of subpaths, each containing a list of fractional k-points.
            For example::

                [ [[0., 0., 0.], [0., 0., 0.5]],
                  [[0.5, 0., 0.], [0.5, 0.5, 0.]] ]

            Will return points along ``0 0 0 -> 0 0 1/2 | 1/2 0 0
            -> 1/2 1/2 0``
        path_labels (:obj:`list`, optional): The k-point labels. These should
            be provided as a :obj:`list` of :obj:`str` for each subpath of the
            overall path. For example::

                [ ['Gamma', 'Z'], ['X', 'M'] ]

            combined with the above example for ``kpt_list`` would indicate the
            path: Gamma -> Z | X -> M. If no labels are provided, letters from
            A -> Z will be used instead.
        eigenvectors (:obj:`bool`, optional): Write the eigenvectors to the
            yaml file.
        dos (str): Path to Phonopy total dos .dat file
        height (:obj:`float`, optional): The height of the plot.
        width (:obj:`float`, optional): The width of the plot.
        ymin (:obj:`float`, optional): The minimum energy on the y-axis.
        ymax (:obj:`float`, optional): The maximum energy on the y-axis.
        style (:obj:`list` or :obj:`str`, optional): (List of) matplotlib style
            specifications, to be composed on top of Sumo base style.
        no_base_style (:obj:`bool`, optional): Prevent use of sumo base style.
            This can make alternative styles behave more predictably.
        image_format (:obj:`str`, optional): The image file format. Can be any
            format supported by matplotlib, including: png, jpg, pdf, and svg.
            Defaults to pdf.
        dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
            the image.
        plt (:obj:`matplotlib.pyplot`, optional): A
            :obj:`matplotlib.pyplot` object to use for plotting.
        fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
            a single font, specified as a :obj:`str`, or several fonts,
            specified as a :obj:`list` of :obj:`str`.

    Returns:
        A matplotlib pyplot object.
    """
    if '.yaml' in filename:
        yaml_file = filename
    elif ('FORCE_SETS' == filename or 'FORCE_CONSTANTS' == filename or
            '.hdf5' in filename):
        try:
            poscar = poscar if poscar else 'POSCAR'
            poscar = Poscar.from_file(poscar)
        except IOError:
            msg = "Cannot find POSCAR file, cannot generate symmetry path."
            logging.error("\n {}".format(msg))
            sys.exit()

        if not dim:
            logging.info("Supercell size (--dim option) not provided.\n"
                         "Attempting to guess supercell dimensions.\n")
            try:
                sposcar = Poscar.from_file("SPOSCAR")
            except IOError:
                msg = "Could not determine supercell size (use --dim flag)."
                logging.error("\n {}".format(msg))
                sys.exit()

            dim = (sposcar.structure.lattice.matrix *
                   poscar.structure.lattice.inv_matrix)

            # round due to numerical noise error
            dim = np.around(dim, 5)

        elif len(dim) == 9:
            dim = np.array(dim).reshape(3, 3)

        elif np.array(dim).shape != (3, 3):
            dim = np.diagflat(dim)

        logging.info("Using supercell with dimensions:")
        logging.info('\t' + str(dim).replace('\n', '\n\t')+'\n')

        factors = {'ev': VaspToEv, 'thz': VaspToTHz, 'mev': VaspToEv * 1000,
                   'cm-1': VaspToCm}

        phonon = load_phonopy(filename, poscar.structure, dim, symprec=symprec,
                              primitive_matrix=primitive_axis,
                              factor=factors[units.lower()],
                              symmetrise=True, born=born,
                              write_fc=False)

        # calculate band structure
        kpath, kpoints, labels = get_path_data(poscar.structure, mode=mode,
                                               symprec=symprec, spg=spg,
                                               kpt_list=kpt_list,
                                               labels=labels, phonopy=True)

        # todo: calculate dos and plot also
        # phonon.set_mesh(mesh, is_gamma_center=False, is_eigenvectors=True,
        #                 is_mesh_symmetry=False)
        # phonon.set_partial_DOS()

        phonon.set_band_structure(
            kpoints, is_eigenvectors=eigenvectors, labels=labels)
        yaml_file = 'sumo_band.yaml'
        phonon._band_structure.write_yaml(filename=yaml_file)

    else:
        msg = "Do not recognise file type of {}".format(filename)
        logging.error("\n {}".format(msg))
        sys.exit()

    save_files = False if plt else True  # don't save if pyplot object provided

    bs = get_ph_bs_symm_line(yaml_file, has_nac=False,
                             labels_dict=kpath.kpoints)

    # Replace dos filename with data array
    if dos is not None:
        if isfile(dos):
            dos = np.genfromtxt(dos, comments='#')
        elif dos:
            phonon.set_mesh(qmesh, is_gamma_center=False, is_eigenvectors=True,
                            is_mesh_symmetry=False)
            phonon.set_total_DOS()
            dos_freq, dos_val = phonon.get_total_DOS()
            dos = np.zeros((len(dos_freq), 2))
            dos[:, 0], dos[:, 1] = dos_freq, dos_val

    plotter = SPhononBSPlotter(bs)
    plt = plotter.get_plot(units=units, ymin=ymin, ymax=ymax, height=height,
                           width=width, plt=plt, fonts=fonts, dos=dos)

    if save_files:
        basename = 'phonon_band.{}'.format(image_format)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename

        if directory:
            filename = os.path.join(directory, filename)

        if dpi is None:
            dpi = rcParams['figure.dpi']
        plt.savefig(filename, format=image_format, dpi=dpi,
                    bbox_inches='tight')

        filename = save_data_files(bs, prefix=prefix, directory=directory)
    else:
        return plt


def save_data_files(bs, prefix=None, directory=None):
    """Write the phonon band structure data files to disk.

    Args:
        bs (:obj:`~pymatgen.phonon.bandstructure.PhononBandStructureSymmLine`):
            The phonon band structure.
        prefix (:obj:`str`, optional): Prefix for data file.
        directory (:obj:`str`, optional): Directory in which to save the data.

    Returns:
        str: The filename of the written data file.
    """
    filename = 'phonon_band.dat'
    filename = '{}_phonon_band.dat'.format(prefix) if prefix else filename
    directory = directory if directory else '.'
    filename = os.path.join(directory, filename)

    with open(filename, 'w') as f:
        header = '#k-distance frequency[THz]\n'
        f.write(header)

        for band in bs.bands:
            for d, e in zip(bs.distance, band):
                f.write('{:.8f} {:.8f}\n'.format(d, e))
            f.write('\n')

    return filename


def _get_parser():
    parser = argparse.ArgumentParser(description="""
    phonon-bandplot is a script to produce publication ready phonon band
    structure diagrams""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filename', default='FORCE_SETS', metavar='F',
                        help="FORCE_SETS, FORCE_CONSTANTS or band.yaml file")
    parser.add_argument('-p', '--prefix', metavar='P',
                        help='prefix for the files generated.')
    parser.add_argument('-d', '--directory', metavar='D',
                        help='output directory for files')
    parser.add_argument('-q', '--qmesh', nargs=3, metavar='N',
                        default=(8, 8, 8),
                        help='q-mesh to use for phonon DOS')
    parser.add_argument('-b', '--born', metavar='B',
                        help='born effective charge file')
    parser.add_argument('-e', '--eigenvectors', action='store_true',
                        help='write the phonon eigenvectors to yaml file')
    parser.add_argument('--dim', nargs='+', metavar='N',
                        help='supercell matrix dimensions')
    parser.add_argument('--poscar', default=None, metavar='POS',
                        help="path to POSCAR file (if FORCE_SETS used)")
    parser.add_argument('--primitive-axis', type=float, nargs=9,
                        metavar='M',
                        dest='primitive_axis', default=None,
                        help='conventional to primitive cell transformation')
    parser.add_argument('--symprec', default=0.01, type=float,
                        help='tolerance for finding symmetry (default: 0.01)')
    parser.add_argument('--units', '-u', metavar='UNITS', default='THz',
                        choices=('THz', 'thz', 'cm-1',
                                 'eV', 'ev', 'meV', 'mev'),
                        help=('choose units of phonon frequency '
                              '(THz, cm-1, eV, meV)'))
    parser.add_argument('--spg', type=str, default=None,
                        help='space group number or symbol')
    parser.add_argument('--density', type=int, default=60,
                        help='k-point density along high-symmetry path')
    parser.add_argument('--seekpath', action='store_true',
                        help='use seekpath to generate the high-symmetry path')
    parser.add_argument('--pymatgen', action='store_true',
                        help='use pymatgen to generate the high-symmetry path')
    parser.add_argument('--cartesian', action='store_true',
                        help='use cartesian k-point coordinates')
    parser.add_argument('--kpoints', type=str, default=None,
                        help=('specify a list of kpoints '
                              '(e.g. "0 0 0, 0.5 0 0")'))
    parser.add_argument('--labels', type=str, default=None,
                        help=('specify the labels for kpoints '
                              r'(e.g. "\Gamma,X")'))
    parser.add_argument('--height', type=float, default=None,
                        help='height of the graph')
    parser.add_argument('--width', type=float, default=None,
                        help='width of the graph')
    parser.add_argument('--ymin', type=float, default=None,
                        help='minimum energy on the y-axis')
    parser.add_argument('--ymax', type=float, default=None,
                        help='maximum energy on the y-axis')
    parser.add_argument('--style', type=str, nargs='+', default=None,
                        help='matplotlib style specifications')
    parser.add_argument('--no-base-style', action='store_true',
                        dest='no_base_style',
                        help='prevent use of sumo base style')
    parser.add_argument('--config', type=str, default=None,
                        help='colour configuration file')
    parser.add_argument('--format', type=str, default='pdf',
                        dest='image_format', metavar='FORMAT',
                        help='image file format (options: pdf, svg, jpg, png)')
    parser.add_argument('--dpi', type=int, default=None,
                        help='pixel density for image file')
    parser.add_argument('--font', default=None, help='font to use')
    parser.add_argument('--dos', nargs='?', type=str,
                        default=None, const=True,
                        help='Phonopy .dat file for phonon DOS')
    return parser


def main():
    args = _get_parser().parse_args()
    logging.basicConfig(filename='sumo-phonon-bandplot.log',
                        level=logging.INFO,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    mode = 'bradcrack'
    if args.seekpath:
        mode = 'seekpath'
    elif args.pymatgen:
        mode = 'pymatgen'

    spg = args.spg
    if args.spg:
        try:
            spg = int(spg)
        except ValueError:
            pass

    kpoints = None
    if args.kpoints:
        kpoints = [[list(map(float, kpt.split())) for kpt in kpts.split(',')]
                   for kpts in args.kpoints.split('|')]
    labels = None
    if args.labels:
        labels = [path.split(',') for path in args.labels.split('|')]

    dim = list(map(float, args.dim)) if args.dim else None

    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")
    warnings.filterwarnings("ignore", category=UnicodeWarning,
                            module="matplotlib")
    warnings.filterwarnings("ignore", category=UserWarning,
                            module="pymatgen")

    if args.primitive_axis:
        pa = np.reshape(args.primitive_axis, (3, 3))
    else:
        pa = None

    phonon_bandplot(args.filename, poscar=args.poscar, prefix=args.prefix,
                    directory=args.directory, dim=dim, born=args.born,
                    qmesh=args.qmesh, primitive_axis=pa, symprec=args.symprec,
                    units=args.units, spg=spg, line_density=args.density,
                    mode=mode, kpt_list=kpoints, labels=labels,
                    height=args.height, width=args.width, ymin=args.ymin,
                    ymax=args.ymax, image_format=args.image_format,
                    style=args.style, no_base_style=args.no_base_style,
                    dpi=args.dpi, fonts=args.font,
                    eigenvectors=args.eigenvectors, dos=args.dos)


if __name__ == "__main__":
    main()
