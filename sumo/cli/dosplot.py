# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
A script to plot density of states diagrams.

TODO:
    * Add ability to scale an orbitals density of states
"""

from __future__ import unicode_literals
from pkg_resources import Requirement, resource_filename

import os
import sys
import logging
import argparse
import warnings
from glob import glob

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import sumo.io.questaal
from sumo.electronic_structure.dos import load_dos, write_files
from sumo.plotting.dos_plotter import SDOSPlotter

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "April 9, 2018"


def dosplot(filename=None, code='vasp', prefix=None, directory=None,
            elements=None, lm_orbitals=None, atoms=None, subplot=False,
            shift=True, total_only=False, plot_total=True, legend_on=True,
            legend_frame_on=False, legend_cutoff=3., gaussian=None, height=6.,
            width=8., xmin=-6., xmax=6., num_columns=2, colours=None, yscale=1,
            xlabel='Energy (eV)', ylabel='Arb. units',
            style=None, no_base_style=False,
            image_format='pdf', dpi=400, plt=None, fonts=None):
    """A script to plot the density of states from a vasprun.xml file.

    Args:
        filename (:obj:`str`, optional): Path to a vasprun.xml file (can be
            gzipped).
        code (:obj:`str`, optional): Electronic structure code used ('vasp' or
            'questaal')
        prefix (:obj:`str`, optional): Prefix for file names.
        directory (:obj:`str`, optional): The directory in which to save files.
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
            the POSCAR. For example, the following will calculate the density
            of states for the first 4 Sn atoms and all O atoms in the
            structure::

                {'Sn': (1, 2, 3, 4), 'O': (, )}

            If ``atoms`` is not set or set to ``None`` then all atomic sites
            for all elements will be considered.
        subplot (:obj:`bool`, optional): Plot the density of states for each
            element on separate subplots. Defaults to ``False``.
        shift (:obj:`bool`, optional): Shift the energies such that the valence
            band maximum (or Fermi level for metals) is at 0 eV. Defaults to
            ``True``.
        total_only (:obj:`bool`, optional): Only extract the total density of
            states. Defaults to ``False``.
        plot_total (:obj:`bool`, optional): Plot the total density of states.
            Defaults to ``True``.
        legend_on (:obj:`bool`, optional): Plot the graph legend. Defaults
            to ``True``.
        legend_frame_on (:obj:`bool`, optional): Plot a frame around the
            graph legend. Defaults to ``False``.
        legend_cutoff (:obj:`float`, optional): The cut-off (in % of the
            maximum density of states within the plotting range) for an
            elemental orbital to be labelled in the legend. This prevents
            the legend from containing labels for orbitals that have very
            little contribution in the plotting range.
        gaussian (:obj:`float`, optional): Broaden the density of states using
            convolution with a gaussian function. This parameter controls the
            sigma or standard deviation of the gaussian distribution.
        height (:obj:`float`, optional): The height of the plot.
        width (:obj:`float`, optional): The width of the plot.
        xmin (:obj:`float`, optional): The minimum energy on the x-axis.
        xmax (:obj:`float`, optional): The maximum energy on the x-axis.
        num_columns (:obj:`int`, optional): The number of columns in the
            legend.
        colours (:obj:`dict`, optional): Use custom colours for specific
            element and orbital combinations. Specified as a :obj:`dict` of
            :obj:`dict` of the colours. For example::

                {
                    'Sn': {'s': 'r', 'p': 'b'},
                    'O': {'s': '#000000'}
                }

            The colour can be a hex code, series of rgb value, or any other
            format supported by matplotlib.
        xlabel (:obj:`str`, optional): Label/units for x-axis (i.e. energy)
        ylabel (:obj:`str`, optional): Label/units for y-axis (i.e. DOS)
        yscale (:obj:`float`, optional): Scaling factor for the y-axis.
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

    if code.lower() == 'vasp':
        if not filename:
            if os.path.exists('vasprun.xml'):
                filename = 'vasprun.xml'
            elif os.path.exists('vasprun.xml.gz'):
                filename = 'vasprun.xml.gz'
            else:
                logging.error('ERROR: No vasprun.xml found!')
                sys.exit()

        dos, pdos = load_dos(filename, elements, lm_orbitals, atoms, gaussian,
                             total_only)
    elif code.lower() == 'questaal':
        if filename:
            pdos_file = filename
            ext = pdos_file.split('.')[-1]
        else:
            pdos_candidates = glob('dos.*')
            for candidate in pdos_candidates:
                if candidate.split('.')[-1] in ('pdf', 'png', 'svg',
                                                'jpg', 'jpeg'):
                    continue
                elif candidate.split('.')[-1].lower() in ('gz', 'z', 'bz2'):
                    pdos_file = candidate
                    ext = candidate.split('.')[-2]
                    break
                else:
                    pdos_file = candidate
                    ext = candidate.split('.')[-1]
                    break
            else:
                raise ValueError("No questaal dos file found")

        if os.path.exists('tdos.{}'.format(ext)):
            tdos_file = 'tdos.{}'.format(ext)
        else:
            tdos_file = None
        if os.path.exists('site.{}'.format(ext)):
            site_file = 'site.{}'.format(ext)
        else:
            site_file = None

        if shift:
            logging.warning("Fermi level shift requested, but not implemented "
                            "for Questaal DOS.")

        dos, pdos = sumo.io.questaal.read_dos(
            pdos_file=pdos_file, tdos_file=tdos_file, site_file=site_file,
            ry=True, gaussian=gaussian, total_only=total_only,
            elements=elements, lm_orbitals=lm_orbitals, atoms=atoms)

    save_files = False if plt else True  # don't save if pyplot object provided

    plotter = SDOSPlotter(dos, pdos)
    plt = plotter.get_plot(subplot=subplot, width=width, height=height,
                           xmin=xmin, xmax=xmax, yscale=yscale,
                           colours=colours, plot_total=plot_total,
                           legend_on=legend_on, num_columns=num_columns,
                           legend_frame_on=legend_frame_on,
                           xlabel=xlabel, ylabel=ylabel,
                           legend_cutoff=legend_cutoff, dpi=dpi, plt=plt,
                           fonts=fonts, style=style,
                           no_base_style=no_base_style)

    if save_files:
        basename = 'dos.{}'.format(image_format)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename
        if directory:
            filename = os.path.join(directory, filename)
        plt.savefig(filename, format=image_format, dpi=dpi,
                    bbox_inches='tight')
        write_files(dos, pdos, prefix=prefix, directory=directory)
    else:
        return plt


def _el_orb(string):
    """Parse the element and orbital argument strings.

    The presence of an element without any orbitals means that we want to plot
    all of its orbitals.

    Args:
        string (str): The element and orbitals as a string, in the form
            ``"C.s.p,O"``.

    Returns:
        dict: The elements and orbitals as a :obj:`dict`. For example::

            {'Bi': ['s', 'px', 'py', 'd']}.

        If an element symbol is included with an empty list, then all orbitals
        for that species are considered.
    """
    el_orbs = {}
    for split in string.split(','):
        orbs = split.split('.')
        orbs = [orbs[0], 's', 'p', 'd', 'f'] if len(orbs) == 1 else orbs
        el_orbs[orbs.pop(0)] = orbs
    return el_orbs


def _atoms(atoms_string):
    """Parse the atom string.

    Args:
        atoms_string (str): The atoms to plot, in the form ``"C.1.2.3,"``.

    Returns:
        dict: The atomic indices over which to sum the DOS. Formatted as::

            {Element: [atom_indices]}.

        Indices are zero indexed for each atomic species. If an element symbol
        is included with an empty list, then all sites for that species are
        considered.
    """
    atoms = {}
    for split in atoms_string.split(','):
        sites = split.split('.')
        el = sites.pop(0)
        sites = list(map(int, sites))
        atoms[el] = np.array(sites) - 1
    return atoms


def _get_parser():
    parser = argparse.ArgumentParser(description="""
    dosplot is a script to produce publication-ready density of
    states diagrams""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filename', help='vasprun.xml file to plot',
                        default=None, metavar='F')
    parser.add_argument('-c', '--code', default='vasp', metavar='C',
                        help=('Input file format: "vasp" (vasprun.xml) or '
                              '"questaal" (opt.ext)'))
    parser.add_argument('-p', '--prefix', metavar='P',
                        help='prefix for the files generated')
    parser.add_argument('-d', '--directory', metavar='D',
                        help='output directory for files')
    parser.add_argument('-e', '--elements', type=_el_orb, metavar='E',
                        help='elemental orbitals to plot (e.g. "C.s.p,O")')
    parser.add_argument('-o', '--orbitals', type=_el_orb, metavar='O',
                        help=('orbitals to split into lm-decomposed '
                              'contributions (e.g. "Ru.d")'))
    parser.add_argument('-a', '--atoms', type=_atoms, metavar='A',
                        help=('atoms to include (e.g. "O.1.2.3,Ru.1.2.3")'))
    parser.add_argument('-s', '--subplot', action='store_true',
                        help='plot each element on separate subplots')
    parser.add_argument('-g', '--gaussian', type=float, metavar='G',
                        help='standard deviation of gaussian broadening')
    parser.add_argument('--columns', type=int, default=2, metavar='N',
                        help='number of columns in the legend')
    parser.add_argument('--legend-cutoff', type=float, default=3,
                        dest='legend_cutoff', metavar='C',
                        help=('cut-off in %% of total DOS that determines if'
                              ' a line is given a label (default: 3)'))
    parser.add_argument('--no-legend', action='store_false', dest='legend',
                        help='hide the plot legend')
    parser.add_argument('--legend-frame', action='store_true',
                        dest='legend_frame',
                        help='display a box around the legend')
    parser.add_argument('--no-shift', action='store_false', dest='shift',
                        help='don\'t shift the VBM/Fermi level to 0 eV')
    parser.add_argument('--total-only', action='store_true', dest='total_only',
                        help='only plot the total density of states')
    parser.add_argument('--no-total', action='store_false', dest='total',
                        help='don\'t plot the total density of states')
    parser.add_argument('--height', type=float, default=None,
                        help='height of the graph')
    parser.add_argument('--width', type=float, default=None,
                        help='width of the graph')
    parser.add_argument('--xmin', type=float, default=-6.,
                        help='minimum energy on the x-axis')
    parser.add_argument('--xmax', type=float, default=6.,
                        help='maximum energy on the x-axis')
    parser.add_argument('--config', type=str, default=None,
                        help='colour configuration file')
    parser.add_argument('--style', type=str, nargs='+', default=None,
                        help='matplotlib style specifications')
    parser.add_argument('--no-base-style', action='store_true',
                        dest='no_base_style',
                        help='prevent use of sumo base style')
    parser.add_argument('--xlabel', type=str, default='Energy (eV)',
                        help='x-axis (i.e. energy) label/units')
    parser.add_argument('--ylabel', type=str, default='Arb. units',
                        help='y-axis (i.e. DOS) label/units')
    parser.add_argument('--yscale', type=float, default=1,
                        help='scaling factor for the y-axis')
    parser.add_argument('--format', type=str, default='pdf',
                        dest='image_format', metavar='FORMAT',
                        help='image file format (options: pdf, svg, jpg, png)')
    parser.add_argument('--dpi', type=int, default=400,
                        help='pixel density for image file')
    parser.add_argument('--font', default=None, help='font to use')
    return parser


def main():
    args = _get_parser().parse_args()
    logging.basicConfig(filename='sumo-dosplot.log', level=logging.INFO,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    if args.config is None:
        config_path = resource_filename(Requirement.parse('sumo'),
                                        'sumo/plotting/orbital_colours.conf')
    else:
        config_path = args.config
    colours = configparser.ConfigParser()
    colours.read(os.path.abspath(config_path))

    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")
    warnings.filterwarnings("ignore", category=UnicodeWarning,
                            module="matplotlib")
    warnings.filterwarnings("ignore", category=UserWarning,
                            module="pymatgen")

    dosplot(filename=args.filename, code=args.code, prefix=args.prefix,
            directory=args.directory, elements=args.elements,
            lm_orbitals=args.orbitals, atoms=args.atoms,
            subplot=args.subplot, shift=args.shift, total_only=args.total_only,
            plot_total=args.total, legend_on=args.legend,
            legend_frame_on=args.legend_frame,
            legend_cutoff=args.legend_cutoff, gaussian=args.gaussian,
            height=args.height, width=args.width, xmin=args.xmin,
            xmax=args.xmax, num_columns=args.columns, colours=colours,
            style=args.style, no_base_style=args.no_base_style,
            xlabel=args.xlabel, ylabel=args.ylabel,
            yscale=args.yscale, image_format=args.image_format, dpi=args.dpi,
            fonts=args.font)


if __name__ == "__main__":
    main()
