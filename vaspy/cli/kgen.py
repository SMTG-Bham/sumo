# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import os
import sys
import math
import errno
import shutil
import string
import logging
import argparse

import numpy as np

from vaspy.electronic_structure.bandstructure import (BradCrackKpath,
            SeekpathKpath, PymatgenKpath, get_kpoints, get_kpoints_from_list)

from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen.symmetry.groups import SpaceGroup

"""
A script to generate KPOINTS files for band structure calculations in VASP
"""

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "July 6, 2017"

# TODO:
#  - Add support for CDML labels
#  - Save Brillouin zone diagram
#  - Return a list of filenames/folders

def kgen(filename='POSCAR', directory=None, make_folders=False, symprec=0.01,
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
        filename (str): The name of the VASP structure file.
        directory (str): The directory in which to output files.
        make_folders (bool): Generate folders and copy in relevant files (INCAR,
            POTCAR, POSCAR, and possibily CHGCAR) from the current directory.
        symprec (float): The precision used for determining cell symmetry.
        kpts_per_split (int): If this is set, the k-points are split into
            seperate k-point files (or folders) each containing the number
            of k-points specified. This is useful for hybrid band structure
            calculations where it is often intractable to calculate all k-points
            in the same calculation.
        ibzkpt (str): If this is set, the generated k-points will be appended
            to this file and given a weight of 0. This is necessary for hybrid
            band structure calculations.
        spg (str or int): The space group international number or symbol to
            override the symmetry determined by spglib. This is not recommended
            and only provided for testing purposes.
        line_density (int): The density of k-points along the path.
        mode (str): Sets the method used to calculate the high-symmetry path.
            Choice of 'bradcrack', 'seekpath', and 'pymatgen'.
        cart_coords (bool): Whether the k-points are provided in cartesian
            or reciprocal coordinates.
        kpt_list (list): Manual list of k-points to use. If kpt_list is set it
            will override the mode selection. Should be formatted as a list of
            subpaths, each containing a list of k-points. For example:
            [[[0., 0., 0.], [0., 0., 0.5]], [[0.5, 0., 0.], [0.5, 0.5, 0.]]]
        labels (list): A list of labels to use along with kpt_list. These should
            be provided as a list of subpaths, each containing a list of labels.
            For example: [['Gamma', 'Z'], ['X', 'M']], combined with the above
            kpt_list would indicate the path: Gamma -> Z | X -> M.
            If no labels are provided, letters from A -> Z will be used instead.
    """
    poscar = Poscar.from_file(filename)
    kpath, kpoints, labels = get_kpt_path(poscar.structure, mode=mode,
                                          symprec=symprec, kpt_list=kpt_list,
                                          labels=labels)

    if not kpt_list and not np.allclose(poscar.structure.lattice.matrix,
                                        kpath.prim.lattice.matrix):
        prim_filename = '{}_prim'.format(os.path.basename(filename))
        kpath.prim.to(filename=prim_filename)

        logging.error("\nWARNING: The input structure does not match the "
                      "expected standard\nprimitive symmetry, the path may be "
                      "incorrect! Use at your own risk.\n\nThe correct symmetry "
                      "primitive structure has been saved as {}.".
                      format(prim_filename))

    ibz = _parse_ibzkpt(ibzkpt)

    if make_folders and ibz and kpts_per_split is None:
        logging.info("\nFound {} total kpoints in path, do you want to "
                     "split them up? (y/n)".format(len(kpoints)))
        if input()[0].lower() == 'y':
            logging.info("How many kpoints per file?")
            kpts_per_split = input()

    write_kpoint_files(filename, kpoints, labels, make_folders=make_folders,
                       ibzkpt=ibz, kpts_per_split=kpts_per_split,
                       directory=directory, cart_coords=cart_coords)

def get_kpt_path(structure, mode='bradcrack', symprec=0.01, spg=None,
                 line_density=60, kpt_list=None, labels=None):
    spg = _get_space_group_object(spg, mode)

    if mode == 'bradcrack':
        kpath = BradCrackKpath(structure, symprec=symprec, spg=spg)
    elif mode == 'seekpath':
        kpath = SeekpathKpath(structure, symprec=symprec)
    elif mode == 'pymatgen':
        kpath = PymatgenKpath(structure, symprec=symprec)

    if kpt_list is not None:
        kpoints, labels, path_str, kpt_dict = get_kpoints_from_list(
            structure, kpt_list, path_labels=labels, line_density=density)
    else:
        kpoints, labels = kpath.get_kpoints(line_density=line_density)
        path_str = kpath.path_string
        kpt_dict = kpath.kpoints

    logging.info('Structure information:'.format(structure.num_sites))
    logging.info('\tSpace group number: {}'.format(kpath._spg_data['number']))

    logging.info('\tInternational symbol: {}'.format(kpath.spg_symbol))
    logging.info('\tLattice type: {}'.format(kpath.lattice_type))

    _print_kpath_information(labels, path_str, kpt_dict)
    return kpath, kpoints, labels

def write_kpoint_files(filename, kpoints, labels, make_folders=False,
                       ibzkpt=None, kpts_per_split=None, directory=None,
                       cart_coords=False):
    """Writes the k-points to KPOINTS files or folders.

    Folders are named as 'split-01', 'split-02', etc ...
    KPOINTS files are named KPOINTS_band_split_01 etc ...

    Args:
        kpoints (list): A list of kpoints.
        labels (list): A list of labels (should have same length as kpoints).
            Label should be set to '' for non-high-symmetry points.
        filename (str): The name of the VASP structure file.
        make_folders (bool): Generate folders and copy in relevant files (INCAR,
            POTCAR, POSCAR, and possibily CHGCAR) from the current directory.
        ibzkpt (Ibzkpt): A pymatgen Ibzkpt object. If this is set, the generated
            k-points will be appended to this file and given a weight of 0.
            This is necessary for hybrid band structure calculations.
        kpts_per_split (int): If this is set, the k-points are split into
            seperate k-point files (or folders) each containing the number
            of k-points specified. This is useful for hybrid band structure
            calculations where it is often intractable to calculate all k-points
            in the same calculation.
        directory (str): The directory in which to output files.
        cart_coords (bool): Whether the k-points are provided in cartesian
            or reciprocal coordinates.
    """
    if kpts_per_split:
        kpt_splits = [kpoints[i:i+kpts_per_split] for
                      i in xrange(0, len(kpoints), kpts_per_split)]
        label_splits = [labels[i:i+kpts_per_split] for
                        i in xrange(0, len(labels), kpts_per_split)]
    else:
        kpt_splits = [kpoints]
        label_splits = [labels]

    if cart_coords:
        coord_type = 'cartesian'
        style = Kpoints.supported_modes.Cartesian
    else:
        coord_type = 'reciprocal'
        style = Kpoints.supported_modes.Reciprocal

    kpt_files = []
    for kpt_split, label_split in zip(kpt_splits, label_splits):
        if ibzkpt:
            # hybrid calculation so set k-point weights to 0
            kpt_weights = ibzkpt.kpts_weights + [0] * len(kpt_split)
            kpt_split = ibzkpt.kpts + kpt_split
            label_split = [''] * len(ibzkpt.labels) + label_split
        else:
            # non-SCF calculation so set k-point weights to 1
            kpt_weights = [1] * len(kpt_split)

        segment = ' -> '.join([label for label in label_split if label])
        kpt_file = Kpoints(comment=segment, num_kpts=len(kpt_split),
                           kpts=kpt_split, kpts_weights=kpt_weights,
                           style=style, coord_type=coord_type,
                           labels=label_split)
        kpt_files.append(kpt_file)

    pad = int(math.floor(math.log10(len(kpt_files)))) + 2
    if make_folders:
        for i, kpt_file in enumerate(kpt_files):
            folder = 'split-{}'.format(str(i+1).zfill(pad))
            if directory:
                folder = os.path.join(directory, folder)

            try:
                os.makedirs(folder)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    logging.error("\nERROR: Folders already exist, won't "
                                  "overwrite.")
                    sys.exit()
                else:
                    raise

            kpt_file.write_file(os.path.join(folder, 'KPOINTS'))
            vasp_files = [filename, "INCAR", "POTCAR", "job"]
            vasp_files += [] if ibzkpt else ['CHGCAR']
            for vasp_file in vasp_files:
                if os.path.isfile(vasp_file):
                    shutil.copyfile(vasp_file, os.path.join(folder, vasp_file))
    else:
        for i, kpt_file in enumerate(kpt_files):
            if len(kpt_files) > 1:
                kpt_filename = 'KPOINTS_band_split_{:0d}'.format(i + 1)
            else:
                kpt_filename = 'KPOINTS_band'
            if directory:
                kpt_filename = os.path.join(directory, kpt_filename)
            kpt_file.write_file(kpt_filename)


def _print_kpath_information(labels, path_str, kpt_dict):
    logging.info('\nk-point path:\n\t{}'.format(path_str))
    logging.info('\nk-points:')
    for label, kpoint in iter(kpt_dict.items()):
        coord_str = ' '.join(['{}'.format(c) for c in kpoint])
        logging.info('\t{}: {}'.format(label, coord_str))
    logging.info('\nk-point label indicies:')
    for i, label in enumerate(labels):
        if label:
            logging.info('\t{}: {}'.format(label, i+1))


def _get_space_group_object(spg, mode):
    if spg and mode != 'bradcrack':
        logging.error("ERROR: Specifying symmetry only supported using "
                      "Bradley and Cracknell path.")
        sys.exit()
    elif spg:
        try:
            if type(spg) is int:
                spg = SpaceGroup.from_int_number(spg)
            else:
                spg = SpaceGroup(spg)
            logging.error("WARNING: Forcing space group not recommended, the"
                          " path is likely\nincorrect. Use at your own risk.\n")
        except ValueError:
            logging.error("ERROR: Space group not recognised.")
            sys.exit()
    return spg


def _parse_ibzkpt(ibzkpt):
    if ibzkpt:
        try:
            ibz = Kpoints.from_file(ibzkpt)
            if ibz.tet_number != 0:
                logging.error('\nERROR: IBZKPT contains tetrahedron information.')
                sys.exit()
        except IOError:
            logging.error('\nERROR: Hybrid specified but no IBZKPT file found.')
            sys.exit()
    else:
        ibz = None
    return ibz


def main():
    parser = argparse.ArgumentParser(description="""
    kgen generates KPOINTS files for running band structure calculations in
    VASP. The high symmetry k-point paths defined in Bradley and Cracknell are
    used by default""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-p', '--poscar', default='POSCAR',
                        help='input VASP structure, default is POSCAR',)
    parser.add_argument('-d', '--directory', type=str, default=None,
                        help='output directory for files')
    parser.add_argument('-f', '--folders', action='store_true',
                        help="""generate calculation folders and copy in
                        necessary files""")
    parser.add_argument('-s', '--split', type=int, default=None,
                        help="number of k-points to include per split")
    parser.add_argument('-y', '--hybrid', default=False, action='store_true',
                        help="""append generated k-points to IBZKPT file with
                        zero weight (needed for hybrid band structures)""")
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

    args = parser.parse_args()
    logging.basicConfig(filename='vaspy-kgen.log', level=logging.DEBUG,
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
        kpoints = [[map(float, kpt.split()) for kpt in kpts.split(',')] for
                   kpts in args.kpoints.split('|')]
    labels = None
    if args.labels:
        labels = [path.split(',') for path in args.labels.split('|')]

    kgen(args.poscar, directory=args.directory, symprec=args.symprec,
         make_folders=args.folders, kpts_per_split=args.split,
         ibzkpt=ibzkpt, spg=spg, density=args.density, mode=mode,
         cart_coords=args.cartesian, kpt_list=kpoints, labels=labels)


if __name__ == "__main__":
    main()
