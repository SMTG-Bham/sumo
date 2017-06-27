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
A script to plot density of states (DOS) diagrams using matplotlib
"""

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "March 13, 2017"


line_width = 1.5
empty_space = 1.05
label_size = 22
col_cycle = colour_cycle()

# TODO:
#   - implement magnify state


def dosplot(filename='vasprun.xml', prefix=None, directory=None, elements=None,
            lm_orbitals=None, atoms=None, subplot=False, shift=True,
            total_only=False, plot_total=True, legend_on=True,
            legend_frame_on=False, legend_cutoff=3., gaussian=None, height=6,
            width=8, xmin=-6, xmax=6, num_columns=2, colours=None, yscale=1,
            image_format='pdf', dpi=400, plot_format='mpl', plt=None):
    """Plot the DOS on a single graph.

    Args:
        filename: A vasprun.xml or DOSCAR file to plot. If DOSCAR selected then
            POSCAR should also be present in same directory.
        prefix: What to prefix all generated files with
        plot_total: Whether or not to plot only the total DOS
        elements: If this is set it gives the specific elements and their
            orbitals to plot
        lm_orbitals: If this is set it gives the specific orbitals to split into
            their components
        atoms: A list of atoms of over which to sum the DOS. The index
            starts at 1. If nothing is specified all the atoms are
            considered. Provided in the form {Element: [atoms]}
        unshifted: Don't shift the DOS to the VBM
        subplot: Plot as subplots rather than a standalone graph
        height: The height of the plot
        width: The width of the plot
        emin: The minimum energy on the x axis
        emax: The maximum energy on the x axis
        ncol: The number of columns of the legend
        tol: The tolerance for finding the VBM
        smear: The amount of Gaussian smearing to apply
        no_total: Whether or not to plot the total DOS
        no_total: Config parser
        colours: A dict of {'Element': colour} where colour is a hex number
    """
    vr = Vasprun(filename)
    band = vr.get_band_structure()
    dos = vr.complete_dos

    if band.is_metal():
        logging.info('System is metallic')
        zero_point = vr.efermi
    else:
        logging.info('Band gap: {:.3f}'.format(band.get_band_gap()['energy']))
        logging.info('DOS band gap: {:.3f}'.format(dos.get_gap()))
        zero_point = band.get_vbm()['energy']

    if shift:
        dos.energies -= zero_point
        if vr.parameters['ISMEAR'] == 0 or vr.parameters['ISMEAR'] == -1:
            dos.energies -= vr.parameters['SIGMA']

    if gaussian:
        dos = dos.get_smeared_vaspdos(gaussian)

    pdos = {}
    if not total_only:
        pdos = get_pdos(dos, lm_orbitals=lm_orbitals, atoms=atoms,
                        elements=elements)

    write_files(dos, pdos, prefix=prefix, directory=directory)

    return plot_figure(dos, pdos, plot_format=plot_format, prefix=prefix,
                       directory=directory, subplot=subplot, width=width,
                       height=height, xmin=xmin, xmax=xmax, yscale=yscale,
                       colours=colours, plot_total=plot_total,
                       legend_on=legend_on, num_columns=num_columns,
                       legend_frame_on=legend_frame_on,
                       legend_cutoff=legend_cutoff,
                       image_format=image_format, dpi=dpi, plt=plt)


def plot_figure(dos, pdos, plot_format='mpl', prefix=None, directory=None,
                subplot=False, width=8, height=6, xmin=-6, xmax=6, yscale=1,
                colours=None, plot_total=True, legend_on=True, num_columns=2,
                legend_frame_on=False, legend_cutoff=3, image_format='pdf',
                dpi=400, plt=None):
    # build a big dictionary of our plotting data then pass to relevant method
    # mask needed to prevent unwanted data in pdf and for finding y limit
    mask = (dos.energies >= xmin - 0.05) & (dos.energies <= xmax + 0.05)
    plot_data = {'mask': mask, 'xmin': xmin, 'xmax': xmax, 'ncol': num_columns,
                 'energies': dos.energies, 'width': width, 'height': height, 
                 'legend_on': legend_on, 'legend_frame_on': legend_frame_on, 
                 'subplot': subplot}
    spins = dos.densities.keys()
    ymax = 0

    if plot_total:
        lines = []
        tdos = {'label': 'Total DOS', 'dens': dos.densities, 'colour': 'k',
                'alpha': 0.15}
        # subplot data formatted as a list of lists of dicts, with each list of
        # dicts being plotted on a seperate graph, if only one list then solo
        # plot
        lines.append([tdos])
        dmax = max([max(d[mask]) for d in dos.densities.values()])
        ymax = dmax if dmax > ymax else ymax
    elif not subplot:
        lines = [[]] # need a blank list to add lines into

    # TODO: Fix broken behaviour if plot_total is off
    cutoff = (legend_cutoff / 100.) * (ymax / 1.05)

    for el, el_pdos in pdos.iteritems():
        el_lines = []
        for orb in sort_orbitals(el_pdos):
            dmax = max([max(d[mask]) for d in el_pdos[orb].densities.values()])
            ymax = dmax if dmax > ymax else ymax
            label = None if dmax < cutoff else '{} ({})'.format(el, orb)
            colour = get_colour_for_element_and_orbital(el, orb, colours)
            el_lines.append({'label': label, 'alpha': 0.25, 'colour': colour,
                             'dens': el_pdos[orb].densities})
        if subplot:
            lines.append(el_lines)
        else:
            lines[0].extend(el_lines)

    ymax = ymax * empty_space / yscale
    ymin = 0 if len(spins) == 1 else -ymax
    plot_data.update({'lines': lines, 'ymax': ymax, 'ymin': ymin})

    if plot_format == 'mpl':
        return _plot_mpl(plot_data, prefix=prefix, directory=directory,
                         image_format=image_format, dpi=dpi, plt=plt)
    elif plot_format == 'xmgrace':
        return _plot_xmgrace(plot_data, prefix=prefix, directory=directory,
                             dpi=dpi)


def _plot_mpl(plot_data, prefix=None, directory=None, image_format='pdf',
              dpi=400, plt=None):

    if plot_data['subplot']:
        nplots = len(plot_data['lines']) + 1
        plt = pretty_subplot(nplots, width=plot_data['width'],
                             height=plot_data['height'], dpi=dpi, plt=plt)
    else:
        plt = pretty_plot(width=plot_data['width'], height=plot_data['height'],
                          dpi=dpi, plt=plt)

    fig = plt.gcf()
    lines = plot_data['lines']
    spins = [Spin.up] if len(lines[0][0]['dens']) == 1 else [Spin.up, Spin.down]
    for i, line_set in enumerate(plot_data['lines']):
        if plot_data['subplot']:
            ax = fig.axes[i]
        else:
            ax = plt.gca()

        for line in line_set:
            for spin in spins:
                if spin == Spin.up:
                    label = line['label']
                    densities = line['dens'][spin]
                elif spin == Spin.down:
                    label = ""
                    densities = -line['dens'][spin]
                ax.fill_between(plot_data['energies'], densities, lw=0,
                                facecolor=line['colour'], alpha=line['alpha'])
                ax.plot(plot_data['energies'], densities, label=label,
                        color=line['colour'], lw=line_width)

        ax.set_ylim(plot_data['ymin'], plot_data['ymax'])
        ax.set_xlim(plot_data['xmin'], plot_data['xmax'])

        ax.tick_params(axis='x', which='both', top='off')
        ax.tick_params(axis='y', which='both', labelleft='off',
                       labelright='off', left='off', right='off')

        loc = 'upper right' if plot_data['subplot'] else 'best'
        ncol = 1 if plot_data['subplot'] else plot_data['ncol']
        if plot_data['legend_on']:
            ax.legend(loc=loc, frameon=plot_data['legend_frame_on'], ncol=ncol,
                      prop={'size': label_size - 3})

    # no add axis labels and sort out ticks
    if plot_data['subplot']:
        fig.text(0.08, 0.5, 'Arb.units', fontsize=label_size, ha='center',
                 va='center', rotation='vertical')
        ax.set_xlabel('Energy (eV)', fontsize=label_size)
        fig.subplots_adjust(hspace=0)
        #plt.tick_params(axis='both', labelsize=label_size)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    else:
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Arb.units')

    plt.tight_layout()
    basename = 'dos.{}'.format(image_format)
    filename = '{}_{}'.format(prefix, basename) if prefix else basename
    if directory:
        filename = os.path.join(directory, filename)

    plt.savefig(filename, format=image_format, dpi=dpi)
    return plt


def _plot_xmgrace(plot_data, prefix=None, directory=None, dpi=400):
    return 0


def el_orb(string):
    """Parse the element and orbital argument strings.

    The presence of an element without any orbitals means that we want to plot
    all of its orbitals.

    Args:
        string: The supplied argument in the form "C.s.p,O"

        Returns:
        The elements and orbitals to plot in the form {Element:[Orbitals]}
    """
    el_orbs = {}
    for split in string.split(','):
        orbs = split.split('.')
        orbs = [orbs[0], 's', 'p', 'd', 'f'] if len(orbs) == 1 else orbs
        el_orbs[orbs.pop(0)] = orbs
    return el_orbs


def get_colour_for_element_and_orbital(element, orbital, colours=None):
    """Choose a colour for a particular elemental orbital.

    If that element is not specified in the colours dictionary, a random colour
    will be generated.

    Args:
        element: The element to choose a colour for.
        orbital: The orbital.
        colours: An dict of {'Element': {'orb': colour}} where colour is a hex
                 number.

    Returns:
        A colour (either the colour name or hex code).
    """
    try:
        return colours.get(element, orbital)
    except (configparser.NoSectionError, configparser.NoOptionError):
        return [a/255. for a in col_cycle.next()]


def atoms(string):
    """Parse the atom string.

    Args:
        string: The supplied argument in the form "C.1.2.3,"

    Returns:
        The elements and orbitals to plot in the form {Element:[Orbitals]}.
        Atoms are zero indexed.
    """
    atoms = {}
    for split in string.split(','):
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
