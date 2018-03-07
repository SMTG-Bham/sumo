# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import os
import sys
import logging
import argparse
import warnings

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

from vaspy.electronic_structure.dos import load_dos, get_pdos, write_files
from vaspy.electronic_structure.plotter import VDOSPlotter

from pymatgen.electronic_structure.core import Spin

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

"""
A script to plot density of states (DOS) diagrams
"""

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "March 13, 2017"


# TODO:
#   - implement magnify state


def dosplot(filename='vasprun.xml', prefix=None, directory=None, elements=None,
            lm_orbitals=None, atoms=None, subplot=False, shift=True,
            total_only=False, plot_total=True, legend_on=True,
            legend_frame_on=False, legend_cutoff=3., gaussian=None, height=6.,
            width=8., xmin=-6., xmax=6., num_columns=2, colours=None, yscale=1,
            image_format='pdf', dpi=400, plt=None, fonts=None):
    """A script to plot the density of states from a vasprun.xml file.

    Args:
        filename (str): A vasprun.xml file to plot (can be gziped).
        prefix (str): A prefix for the files generated.
        directory (str): Specify a directory in which the files are saved.
        elements (dict): A dict of element names specifying which orbitals to
            plot. For example {'Bi': ['s', 'px', 'py', 'd']}. If an element
            symbol is included with an empty list, then all orbitals for that
            species are considered. If set to None then all orbitals for all
            elements are considered.
        lm_orbitals (dict): A list of orbitals for which the lm decomposed
            contributions should be calculated, in the form {Element: [orbs]}
        atoms (dict): A dictionary containing a list of atomic indicies over
            which to sum the DOS, provided as {Element: [atom_indicies]}.
            Indicies are zero indexed for each atomic species. If an element
            symbol is included with an empty list, then all sites for that
            species are considered. If set to None then all sites for all
            elements are considered.
        subplot (bool): Split the plot up into separate plots for each element.
        shift (bool): Shift the energies such that the valence band maximum is
            at 0 eV.
        total_only (bool): Only plot the total density of states.
        plot_total (bool): Whether or not to plot total DOS.
        legend_on (bool): Whether or not to plot the graph legend.
        legend_frame_on (bool): Whether or not to plot the graph legend frame.
        legend_cutoff (int): The cut-off (in % of maximum DOS plotted) for a
            elemental/orbital DOS label to appear in the legend.
        gaussian (float): The sigma of the Gaussian broadening to apply (usually
            controlled by the SIGMA flag in VASP).
        height (float): The height of the graph (matplotlib only).
        width (float): The width of the graph (matplotlib only).
        xmin (float): The minimum energy to plot.
        xmax (float): The maximum energy to plot.
        num_columns (int): The number of columns in the legend.
        colours (dict): Specify custom colours as {'Element': colour} where
            colour is a hex number.
        yscale (dict): Scaling factor for the y-axis.
        image_format (str): The image file format (matplotlib only). Can be
            any format supported by matplot, including: png, jpg, pdf, and svg.
        dpi (int): The dots-per-inch (pixel density) for the image.
        plt (pyplot object): Matplotlib pyplot object to use for plotting.
            If plt is set then no files will be written.
        fonts (list): List of fonts to use.

    Returns:
        A matplotlib pyplot object.
    """
    dos, pdos = load_dos(filename, elements, lm_orbitals, atoms, gaussian,
                         total_only)

    save_files = False if plt else True  # don't save if pyplot object provided

    plotter = VDOSPlotter(dos, pdos)
    plt = plotter.get_plot(subplot=subplot, width=width, height=height,
                           xmin=xmin, xmax=xmax, yscale=yscale, colours=colours,
                           plot_total=plot_total, legend_on=legend_on,
                           num_columns=num_columns,
                           legend_frame_on=legend_frame_on,
                           legend_cutoff=legend_cutoff, dpi=dpi, plt=plt,
                           fonts=fonts)

    if save_files:
        basename = 'dos.{}'.format(image_format)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename
        if directory:
            filename = os.path.join(directory, filename)
        plt.savefig(filename, format=image_format, dpi=dpi)
        write_files(dos, pdos, prefix=prefix, directory=directory)
    else:
        return plt


def el_orb(string):
    """Parse the element and orbital argument strings.

    The presence of an element without any orbitals means that we want to plot
    all of its orbitals.

    Args:
        string (str): The supplied argument in the form "C.s.p,O".

    Returns:
        A dict of element names specifying which orbitals to plot. For example
        {'Bi': ['s', 'px', 'py', 'd']}. If an element symbol is included with
        an empty list, then all orbitals for that species are considered.
    """
    el_orbs = {}
    for split in string.split(','):
        orbs = split.split('.')
        orbs = [orbs[0], 's', 'p', 'd', 'f'] if len(orbs) == 1 else orbs
        el_orbs[orbs.pop(0)] = orbs
    return el_orbs


def atoms(atoms_string):
    """Parse the atom string.

    Args:
        atoms_string (str): The supplied argument in the form "C.1.2.3,".

    Returns:
        A dictionary containing a list of atomic indicies over which to sum
        the DOS, provided as {Element: [atom_indicies]}. Indicies are zero
        indexed for each atomic species. If an element symbol is included with
        an empty list, then all sites for that species are considered.
    """
    atoms = {}
    for split in atoms_string.split(','):
        sites = split.split('.')
        el = sites.pop(0)
        sites = map(int, sites)
        atoms[el] = np.array(sites) - 1
    return atoms


def main():
    parser = argparse.ArgumentParser(description="""
    dosplot is a convenient script to help make publication ready density of
    states diagrams.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filename', help='vasprun.xml file to plot',
                        default='vasprun.xml')
    parser.add_argument('-p', '--prefix', help='prefix for the files generated')
    parser.add_argument('-d', '--directory', help='output directory for files')
    parser.add_argument('-e', '--elements', type=el_orb, help="""Choose the
                        elements to plot. These should be listed using the
                        symbols from the POSCAR and seperated via commas.
                        Specific orbitals can be chosen by adding the orbitals
                        after the element by using a period as the seperator.
                        For example, to plot the carbon s and p, and all the
                        oxygen orbitals, the command would be "-e C.s.p,O".""")
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
    parser.add_argument('--height', type=float, default=6.,
                        help='The height of the graph')
    parser.add_argument('--width', type=float, default=8.,
                        help='The width of the graph')
    parser.add_argument('--xmin', type=float, default=-6.,
                        help='The minimum energy on the x axis')
    parser.add_argument('--xmax', type=float, default=6.,
                        help='The maximum energy on the x axis')
    parser.add_argument('-c', '--columns', type=int, default=2,
                        help='The number of columns in the legend')
    parser.add_argument('--config', type=str, default=None,
                        help='Colour configuration file')
    parser.add_argument('--yscale', type=float, default=1,
                        help='Scaling factor for the y axis')
    parser.add_argument('--format', type=str, default='pdf',
                        dest='image_format',
                        help='select image format from pdf, svg, jpg, & png')
    parser.add_argument('--dpi', type=int, default=400,
                        help='pixel density for generated images')
    parser.add_argument('--font', default=None, help='Font to use.')

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
            fonts=[args.font])


if __name__ == "__main__":
    main()
