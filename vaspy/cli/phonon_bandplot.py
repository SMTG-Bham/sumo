# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import os
import sys
import logging
import argparse
import numpy as np

import warnings
warnings.filterwarnings("ignore", category=FutureWarning,
                        module="h5py")
import h5py

import matplotlib as mpl
mpl.use('Agg')

from phonopy.units import VaspToTHz

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.phonopy import get_ph_bs_symm_line

from vaspy.phonon.phonopy import load_phonopy
from vaspy.plotting.phonon_bs_plotter import VPhononBSPlotter
from vaspy.symmetry import BradCrackKpath, SeekpathKpath, PymatgenKpath
from vaspy.symmetry.kpoints import get_kpoints, get_path_data


"""
A script to plot phonon band structure diagrams
"""


__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Jan 17, 2018"

# TODO:
# - automatically plot dos if present in band.yaml
# - primitive axis determination (try symmetrise->spglib->PA then apply
#   transform on original cell and see if it works.
# - make band structure from vasprun displacement/DFPT files
# - deal with magnetic moments
# - Read FORCE_CONSANTS or force_constants.hdf5
# - change frequency unit
# - read settings from phonopy config file

def phonon_bandplot(filename, poscar=None, prefix=None, directory=None, dim=None,
                    born=None, qmesh=None, spg=None, line_density=60,
                    symprec=0.01, mode='bradcrack', kpt_list=None,
                    eigenvectors=False,
                    labels=None, height=6., width=6., ymin=None, ymax=None,
                    image_format='pdf', dpi=400, plt=None, fonts=None):
    """A script to plot phonon band structure diagrams.

    Args:
        filename (str): A vasprun.xml file to plot (can be gziped).
        prefix (str): A prefix for the files generated.
        directory (str): Specify a directory in which the files are saved.
        born (str):
        qmesh (list):
        spg (str or int): The space group international number or symbol to
            override the symmetry determined by spglib. This is not recommended
            and only provided for testing purposes.
        line_density (int): The density of k-points along the path.
        mode (str): Sets the method used to calculate the high-symmetry path.
            Choice of 'bradcrack', 'seekpath', and 'pymatgen'.
        kpt_list (list): Manual list of k-points to use. If kpt_list is set it
            will override the mode selection. Should be formatted as a list of
            subpaths, each containing a list of k-points. For example:
            [[[0., 0., 0.], [0., 0., 0.5]], [[0.5, 0., 0.], [0.5, 0.5, 0.]]]
        labels (list): A list of labels to use along with kpt_list. These should
            be provided as a list of subpaths, each containing a list of labels.
            For example: [['Gamma', 'Z'], ['X', 'M']], combined with the above
            kpt_list would indicate the path: Gamma -> Z | X -> M.
            If no labels are provided, letters from A -> Z will be used instead.
        eigenvectors (:obj:`bool`, optional): Write the eigenvectors to the
            yaml file.
        height (float): The height of the graph (matplotlib only).
        width (float): The width of the graph (matplotlib only).
        ymin (float): The minimum energy to plot.
        ymax (float): The maximum energy to plot.
        image_format (str): The image file format (matplotlib only). Can be
            any format supported by matplot, including: png, jpg, pdf, and svg.
        dpi (int): The dots-per-inch (pixel density) for the image.
        plt (pyplot object): Matplotlib pyplot object to use for plotting.
            If plt is set then no files will be written.
        fonts (list): List of fonts to use.

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
            try:
                sposcar = Poscar.from_file("SPOSCAR")
            except IOError:
                msg = "Could not determine supercell size (use --dim flag)."
                logging.error("\n {}".format(msg))
                sys.exit()

            dim = sposcar.structure.lattice.matrix * poscar.structure.lattice.inv_matrix
            # round due to numerical noise error
            dim = np.around(dim, 5)

        elif np.array(dim).shape != (3, 3):
            dim = np.diagflat(dim)

        # print dim to user
        phonon = load_phonopy(filename, poscar.structure, dim, symprec=symprec,
                              primitive_matrix=None, factor=VaspToTHz, symmetrise=True,
                              born=born, write_fc=False)

        # calculate band structure
        kpath, kpoints, labels = get_path_data(poscar.structure, mode=mode,
                                           symprec=symprec, kpt_list=kpt_list,
                                           labels=labels, phonopy=True)

        #phonon.set_mesh(mesh, is_gamma_center=False, is_eigenvectors=True,
        #                is_mesh_symmetry=False)
        #phonon.set_partial_DOS()
        phonon.set_band_structure(kpoints, is_eigenvectors=eigenvectors)
        yaml_file = 'vaspy_band.yaml'
        phonon._band_structure.write_yaml(labels=labels, filename=yaml_file)

    else:
        msg = "Do not recognise file type of {}".format(filename)
        logging.error("\n {}".format(msg))
        sys.exit()

    save_files = False if plt else True  # don't save if pyplot object provided

    bs = get_ph_bs_symm_line(yaml_file, has_nac=False,
                             labels_dict=kpath.kpoints)
    plotter = VPhononBSPlotter(bs)
    plt = plotter.get_plot(ymin=ymin, ymax=ymax, height=height, width=width,
                           plt=plt, dos_plotter=None, dos_options=None,
                           fonts=fonts)
    if save_files:
        basename = 'phonon_band.{}'.format(image_format)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename
        if directory:
            filename = os.path.join(directory, filename)
        plt.savefig(filename, format=image_format, dpi=dpi, bbox_inches='tight')

        filename='{}_phonon_band.dat'.format(prefix) if prefix else 'phonon_band.dat'
        with open(filename, 'w') as f:
            header = '#k-distance frequency[THz]\n'
            f.write(header)
            for band in bs.bands:
                for d, e in zip(bs.distance, band):
                    f.write('{:.8f} {:.8f}\n'.format(d, e))
                f.write('\n')
    else:
        return plt

def main():
    parser = argparse.ArgumentParser(description="""
    dosplot is a convenient script to help make publication ready density of
    states diagrams.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filename', default='FORCE_SETS',
                        help="phonopy band.conf or FORCE_SETS file")
    parser.add_argument('--poscar', default=None,
                        help="Path to POSCAR file, needed if FORCE_SETS used.")
    parser.add_argument('-p', '--prefix', help='Prefix for the files generated.')
    parser.add_argument('-d', '--directory', help='output directory for files')
    parser.add_argument('--dim', nargs='+',
                        help='Supercell matrix dimensions')
    parser.add_argument('-q', '--qmesh', nargs=3,
                        help='q-mesh to use for phonon DOS')
    parser.add_argument('-b', '--born', nargs='*',
                        help="""File containing born effective charges. Can
                        be the output from phonopy-vasp-born, or a vasprun.xml
                        file""")
    parser.add_argument('--symprec', default=0.01, type=float,
                        help='tolerance for finding symmetry, default is 0.01')
    parser.add_argument('--spg', type=str, default=None,
                        help='space group number to override detected symmetry')
    parser.add_argument('--density', type=int, default=60,
                        help='k-point density along high symmetry lines')
    parser.add_argument('--seekpath', action='store_true',
                        help='use seekpath to generate the high-symmetry path')
    parser.add_argument('--pymatgen', action='store_true',
                        help='use pymatgen to generate the high-symmetry path')
    parser.add_argument('--cartesian', action='store_true',
                        help='use cartesian rather than fractional coordinates')
    parser.add_argument('--kpoints', type=str, default=None,
                        help="""specify a list of kpoints manually, written as
                        --kpoints '0 0 0, 0.5 0 0'""")
    parser.add_argument('--labels', type=str, default=None,
                        help="""specify the labels for manual kpoints, written
                        as --labels '\Gamma,X'""")
    parser.add_argument('-e', '--eigenvectors', action='store_true',
                        help='Write the phonon eigenvectors to the yaml file.')
    parser.add_argument('--height', type=float, default=6.0,
                        help='The height of the graph')
    parser.add_argument('--width', type=float, default=6.0,
                        help='The width of the graph')
    parser.add_argument('--ymin', type=float, default=None,
                        help='The minimum energy on the y axis')
    parser.add_argument('--ymax', type=float, default=None,
                        help='The maximum energy on the y axis')
    parser.add_argument('--format', type=str, default='pdf',
                        dest='image_format',
                        help='select image format from pdf, svg, jpg, & png')
    parser.add_argument('--dpi', type=int, default=400,
                        help='pixel density for generated images')
    parser.add_argument('--font', default=None, help='Font to use.')

    args = parser.parse_args()
    logging.basicConfig(filename='vaspy-phonon-bandplot.log',
                        level=logging.DEBUG,
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
        kpoints = [[list(map(float, kpt.split())) for kpt in kpts.split(',')] for
                   kpts in args.kpoints.split('|')]
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

    phonon_bandplot(args.filename, poscar=args.poscar, prefix=args.prefix, directory=args.directory, dim=dim,
                    born=args.born, qmesh=args.qmesh, spg=spg, line_density=args.density, mode=mode, kpt_list=kpoints,
             labels=labels, height=args.height, width=args.width, ymin=args.ymin, ymax=args.ymax,
             image_format=args.image_format, dpi=args.dpi, fonts=[args.font],
                    eigenvectors=args.eigenvectors)

if __name__ == "__main__":
    main()
