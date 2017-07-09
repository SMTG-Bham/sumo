#!/usr/bin/env python
# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import os
import sys
import logging
import argparse
import itertools
import warnings

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

from vaspy.electronic_structure.dos import sort_orbitals, get_pdos, write_files
from vaspy.misc.plotting import pretty_plot, pretty_subplot, colour_cycle

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

"""
A script to generate KPOINTS files for band structure calculations in VASP
"""

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "July 6, 2017"

# load structure file

def kgen(filename='vasprun.xml', symprec=1e-3, mode='seekpath', hybrid=False,
        kpts_per_folder=None, make_folders=False, spg=None):




def main():
    parser = argparse.ArgumentParser(description="""
    dosplot is a convenient script to help make publication ready density of
    states diagrams.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filename', default='POSCAR'
                        help='input VASP structure, default is POSCAR',)
    parser.add_argument('-t', '--tolerance', default=1e-3, type=float,
                        help='tolerance for finding symmetry')
    parser.add_argument('-o', '--output', help='output directory for files')
    parser.add_argument('-f', '--folders', action='store_true',
                        help="""generate calculation folders and copy relevant 
                        files""")
    parser.add_argument('-o', '--orbitals', type=el_orb, help="""Choose the
                        orbital to split. This should be listed as the element
                        (using the symbol from the POSCAR) and the orbital
                        seperated by a period. For example to plot the oxygen
                        split d orbitals, the command would be "-o O.d". More
                        than one split orbital and element can be added using
                        the notation described for adding more elements.""")
    parser.add_argument('-a', '--atoms', type=atoms, help="""Choose which atoms
                        to calculate the DOS for. This should be listed as the
                        element (using the symbol from the POSCAR) and the atoms
                        seperated by a period. For example to plot the oxygen 1,
                        2 and 3 atoms, the command would be "-a O.1.2.3". The
                        atom indicies start at 1 (as in the VASP output). You
                        can specify a range to avoid typing all the numbers
                        out, e.g. the previous command can be written "-a
                        O.1-3". To select all the atoms of an element just
                        include the element symbol with no numbers after it,
                         e.g. "-a Ru" will include all the Ru atoms. If
                        an element is not specified then it will not be
                        included in the DOS. More than one element can be added
                        using the notation described above for adding more
                        elements.""")
    parser.add_argument('-s', '--subplot', action='store_true', help="""Plot the
                        DOS as a series of subplots rather than on a single
                        graph. The height and width arguments for this program
                        refer to the dimensions of a single subplot rather than
                        the overall figure.""")
    parser.add_argument('--no-shift', action='store_false', dest='shift',
                        help='Don\'t shift the graph so that the VBM is at 0')
    parser.add_argument('--total-only', action='store_true', dest='total_only',
                        help='Only plot the total DOS')
    parser.add_argument('--no-total', action='store_false', dest='total',
                        help='Don\'t plot the total DOS')
    parser.add_argument('--no-legend', action='store_false', dest='legend',
                        help='Don\'t display the plot legend')
    parser.add_argument('--legend-frame', action='store_true',
                        dest='legend_frame',
                        help='Display a frame box around the graph legend.')
    parser.add_argument('--legend-cutoff', type=float, default=3,
                        dest='legend_cutoff',
                        help="""Cut-off in %% of total DOS in visible range that
                        determines if a line is given a label. Set to 0 to label
                        all lines. Default is 3 %%""")
    parser.add_argument('-g', '--gaussian', type=float,
                        help='Amount of gaussian broadening to apply')
    parser.add_argument('--height', type=float, default=6,
                        help='The height of the graph')
    parser.add_argument('--width', type=float, default=8,
                        help='The width of the graph')
    parser.add_argument('--xmin', type=float, default=-6.0,
                        help='The minimum energy on the x axis')
    parser.add_argument('--xmax', type=float, default=6.0,
                        help='The maximum energy on the x axis')
    parser.add_argument('-c', '--columns', type=int, default=2,
                        help='The number of columns in the legend')
    parser.add_argument('--config', type=str, default=None,
                        help='Colour configuration file')
    parser.add_argument('--yscale', type=float, default=1,
                        help='Scaling factor for the y axis')
    parser.add_argument('--xmgrace', action='store_true',
                        help='plot using xmgrace instead of matplotlib')
    parser.add_argument('--image_format', type=str, default='pdf',
                        help='select image format from pdf, svg, jpg, & png')
    parser.add_argument('--dpi', type=int,
                        help='pixel density for generated images')

    args = parser.parse_args()
    logging.basicConfig(filename='vaspy-dosplot.log', level=logging.DEBUG,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    if args.config is None:
        config_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'default_colours.ini')
    else:
        config_path = args.config
    colours = configparser.ConfigParser()
    colours.read(os.path.abspath(config_path))

    if args.xmgrace:
        plot_format = 'xmgrace'
    else:
        plot_format = 'mpl'
        warnings.filterwarnings("ignore", category=UserWarning, 
                                module="matplotlib")

    dosplot(filename=args.filename, prefix=args.prefix, directory=args.directory,
            elements=args.elements, lm_orbitals=args.orbitals, atoms=args.atoms,
            subplot=args.subplot, shift=args.shift, total_only=args.total_only, 
            plot_total=args.total, legend_on=args.legend,
            legend_frame_on=args.legend_frame,
            legend_cutoff=args.legend_cutoff, gaussian=args.gaussian,
            height=args.height, width=args.width, xmin=args.xmin,
            xmax=args.xmax, num_columns=args.columns, colours=colours,
            yscale=args.yscale, image_format=args.image_format, dpi=args.dpi,
            plot_format=plot_format)


if __name__ == "__main__":
    main()
