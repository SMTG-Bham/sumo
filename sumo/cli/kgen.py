# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
A script to generate KPOINTS files for band structure calculations in VASP

TODO:
  * Add support for CDML labels
  * Save Brillouin zone diagram
  * Return a list of filenames/folders
"""

from __future__ import unicode_literals

import os
import sys
import logging
import argparse

import numpy as np

from sumo.symmetry.kpoints import get_path_data, write_kpoint_files
from pymatgen.io.vasp.inputs import Poscar, Kpoints
import sumo.io.questaal
from sumo.io.questaal import QuestaalInit, QuestaalSite


__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "July 6, 2017"


def kgen(filename='POSCAR', code='vasp',
         directory=None, make_folders=False, symprec=0.01,
         kpts_per_split=None, ibzkpt=None, spg=None, density=60,
         mode='bradcrack', cart_coords=False, kpt_list=None, labels=None):
    """Generate KPOINTS files for VASP band structure calculations.

    This script provides a wrapper around several frameworks used to generate
    k-points along a high-symmetry path. The paths found in Bradley and
    Cracknell, SeeK-path, and pymatgen are all supported.

    It is important to note that the standard primitive cell symmetry is
    different between SeeK-path and pymatgen. If the correct the structure
    is not used, the high-symmetry points (and band path) may be invalid.

    Args:
        filename (:obj:`str`, optional): Path to VASP structure file. Default
            is ``POSCAR``.
        code (:obj:`str`, optional): Calculation type. Default is 'vasp';
            'questaal' also supported.
        directory (:obj:`str`, optional): The output file directory.
        make_folders (:obj:`bool`, optional): Generate folders and copy in
            required files (INCAR, POTCAR, POSCAR, and possibly CHGCAR) from
            the current directory.
        symprec (:obj:`float`, optional): The precision used for determining
            the cell symmetry.
        kpts_per_split (:obj:`int`, optional): If set, the k-points are split
            into separate k-point files (or folders) each containing the number
            of k-points specified. This is useful for hybrid band structure
            calculations where it is often intractable to calculate all
            k-points in the same calculation.
        ibzkpt (:obj:`str`, optional): Path to IBZKPT file. If set, the
            generated k-points will be appended to the k-points in this file
            and given a weight of 0. This is necessary for hybrid band
            structure calculations.
        spg (:obj:`str` or :obj:`int`, optional): The space group international
            number or symbol to override the symmetry determined by spglib.
            This is not recommended and only provided for testing purposes.
            This option will only take effect when ``mode = 'bradcrack'``.
        line_density (:obj:`int`, optional): Density of k-points along the
            path.
        mode (:obj:`str`, optional): Method used for calculating the
            high-symmetry path. The options are:

            bradcrack
                Use the paths from Bradley and Cracknell. See [brad]_.

            pymatgen
                Use the paths from pymatgen. See [curt]_.

            seekpath
                Use the paths from SeeK-path. See [seek]_.

        cart_coords (:obj:`bool`, optional): Whether the k-points are returned
            in cartesian or reciprocal coordinates. Defaults to ``False``
            (fractional coordinates).
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
            A -> Z will be used instead. If a label begins with '@' it will be
            concealed when plotting with sumo-bandplot.
    """
    if code.lower() == 'vasp':
        structure = Poscar.from_file(filename).structure
    elif code.lower() == 'questaal':
        if filename.split('.')[0] == 'site':
            site = QuestaalSite.from_file(filename)
            structure = site.structure
            alat = site.alat
        else:
            structure = QuestaalInit.from_file(filename).structure
            alat = None
    else:
        raise ValueError('Code "{0}" not recognized.'.format(code))

    kpath, kpoints, labels = get_path_data(structure, mode=mode,
                                           symprec=symprec, kpt_list=kpt_list,
                                           labels=labels, spg=spg,
                                           line_density=density,
                                           cart_coords=cart_coords)

    logging.info('\nk-point label indices:')
    for i, label in enumerate(labels):
        if label:
            logging.info('\t{}: {}'.format(label, i+1))

    if not kpt_list and not np.allclose(structure.lattice.matrix,
                                        kpath.prim.lattice.matrix):
        prim_filename = '{}_prim'.format(os.path.basename(filename))
        if code.lower() == 'questaal':
            QuestaalInit.from_structure(kpath.prim).to_file(prim_filename)
        else:
            kpath.prim.to(filename=prim_filename)

        logging.error("\nWARNING: The input structure does not match the "
                      "expected standard\nprimitive symmetry, the path may be "
                      "incorrect! Use at your own risk.\n\nThe correct "
                      "symmetry primitive structure has been saved as {}.".
                      format(prim_filename))

    ibz = _parse_ibzkpt(ibzkpt)

    if make_folders and ibz and kpts_per_split is None:
        logging.info("\nFound {} total kpoints in path, do you want to "
                     "split them up? (y/n)".format(len(kpoints)))
        if input()[0].lower() == 'y':
            logging.info("How many kpoints per file?")
            kpts_per_split = int(input())

    if code.lower() == 'vasp':
        write_kpoint_files(filename, kpoints, labels, make_folders=make_folders,
                           ibzkpt=ibz, kpts_per_split=kpts_per_split,
                           directory=directory, cart_coords=cart_coords)

    elif code.lower() == 'questaal':
        if cart_coords:
            kpoints = [kpoint / (2 * np.pi) for kpoint in kpoints]
            if alat is not None:
                logging.info("Multiplying kpoint values by ALAT = {} Bohr".format(alat))
                _bohr_to_angstrom = 0.5291772
                kpoints = [kpoint * alat * _bohr_to_angstrom for kpoint in kpoints]
        sumo.io.questaal.write_kpoint_files(filename, kpoints, labels,
                                            make_folders=make_folders,
                                            directory=directory,
                                            cart_coords=cart_coords)

def _parse_ibzkpt(ibzkpt):
    if ibzkpt:
        try:
            ibz = Kpoints.from_file(ibzkpt)
            if ibz.tet_number != 0:
                logging.error('\nERROR: IBZKPT contains tetrahedron '
                              'information.')
                sys.exit()
        except IOError:
            logging.error('\nERROR: Hybrid specified but no IBZKPT file '
                          'found.')
            sys.exit()
    else:
        ibz = None
    return ibz


def _get_parser():
    parser = argparse.ArgumentParser(description="""
    kgen generates KPOINTS files for running band structure calculations in
    VASP""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-p', '--poscar', default='POSCAR', metavar='P',
                        help='input structure file (default: POSCAR)',)
    parser.add_argument('-c', '--code', default='vasp',
                        help='Electronic structure code (default: vasp).'
                             '"questaal" also supported.')
    parser.add_argument('-d', '--directory', type=str, default=None,
                        metavar='D', help='output directory for files')
    parser.add_argument('-s', '--split', type=int, default=None, metavar='N',
                        help='number of k-points to include per file')
    parser.add_argument('-f', '--folders', action='store_true',
                        help='create folders and copy in necessary files')
    parser.add_argument('-y', '--hybrid', default=False, action='store_true',
                        help='append k-points to IBZKPT file with zero weight')
    parser.add_argument('--symprec', default=0.01, type=float,
                        help='tolerance for finding symmetry (default: 0.01)')
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
    return parser


def main():
    args = _get_parser().parse_args()

    logging.basicConfig(filename='sumo-kgen.log', level=logging.DEBUG,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    mode = 'bradcrack'
    if args.seekpath:
        mode = 'seekpath'
    elif args.pymatgen:
        mode = 'pymatgen'

    ibzkpt = 'IBZKPT' if args.hybrid else None

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

    kgen(args.poscar, code=args.code,
         directory=args.directory, symprec=args.symprec,
         make_folders=args.folders, kpts_per_split=args.split,
         ibzkpt=ibzkpt, spg=spg, density=args.density, mode=mode,
         cart_coords=args.cartesian, kpt_list=kpoints, labels=labels)


if __name__ == "__main__":
    main()
