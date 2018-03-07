# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import os
import sys
import glob
import logging
import argparse
import warnings

import numpy as np
import matplotlib as mpl
mpl.use('Agg')

from pkg_resources import Requirement, resource_filename

from vaspy.plotting.bs_plotter import VBSPlotter
from vaspy.plotting.dos_plotter import VDOSPlotter
from vaspy.electronic_structure.dos import load_dos, get_pdos
from vaspy.cli.dosplot import atoms, el_orb

from pymatgen.io.vasp.outputs import BSVasprun
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import \
    get_reconstructed_band_structure

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

"""
A script to plot band structure diagrams
"""

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "July 18, 2017"


line_width = 1.5
empty_space = 1.05
label_size = 22


def bandplot(filenames=None, prefix=None, directory=None, vbm_cbm_marker=False,
             project=None, project_rgb=None, dos_file=None, elements=None,
             lm_orbitals=None, atoms=None, total_only=False,
             plot_total=True, legend_cutoff=3, gaussian=None, height=6.,
             width=6., ymin=-6., ymax=6., colours=None, yscale=1,
             image_format='pdf', dpi=400, plt=None, fonts=None):
    if not filenames:
        folders = glob.glob('split-*')
        folders = sorted(folders) if folders else ['.']
        filenames = []
        for fol in folders:
            vr_file = os.path.join(fol, 'vasprun.xml')
            if os.path.exists(vr_file):
                filenames.append(vr_file)
            else:
                logging.error('ERROR: No vasprun.xml found in {}!'.format(fol))
                sys.exit()

    parse_projected = True if project else False
    bandstructures = []
    for vr_file in filenames:
        vr = BSVasprun(vr_file, parse_projected_eigen=parse_projected)
        bs = vr.get_band_structure(line_mode=True)
        bandstructures.append(bs)
    bs = get_reconstructed_band_structure(bandstructures)

    #if project and dos_file:
    #    logging.error('ERROR: Plotting projected band structure with DOS not '
    #                  'supported.\nPlease use --projected-rgb option.')
    #    sys.exit()

    save_files = False if plt else True  # don't save if pyplot object provided

    dos_plotter = None
    dos_opts = None
    if dos_file:
        dos, pdos = load_dos(dos_file, elements, lm_orbitals, atoms, gaussian,
                             total_only)
        dos_plotter = VDOSPlotter(dos, pdos)
        dos_opts = {'plot_total': plot_total, 'legend_cutoff': legend_cutoff,
                    'colours': colours, 'yscale': yscale}

    plotter = VBSPlotter(bs)
    if project:
        elemental_orbitals = [('O', 'p'), ('Bi', 'p'), ('I', 'p')]
        plt = plotter.get_projected_rgb_plot(elemental_orbitals, zero_to_efermi=True,
                              ymin=ymin, ymax=ymax, height=height, width=width,
                              vbm_cbm_marker=vbm_cbm_marker, plt=plt,
                              dos_plotter=dos_plotter, dos_options=dos_opts)
    elif project:
        raise NotImplementedError('projected band structure plotting not yet '
                                  'supported')
    else:
        plt = plotter.get_plot(zero_to_efermi=True, ymin=ymin, ymax=ymax,
                               height=height, width=width,
                               vbm_cbm_marker=vbm_cbm_marker, plt=plt,
                               dos_plotter=dos_plotter, dos_options=dos_opts,
                               fonts=fonts)

    if save_files:
        basename = 'band.{}'.format(image_format)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename
        if directory:
            filename = os.path.join(directory, filename)
        plt.savefig(filename, format=image_format, dpi=dpi, bbox_inches='tight')

        # TODO: save bandstructure dat file properly (spin polarized case
        # and use numpy). This will currently append to the file, even if
        # it already exists.
        filename='{}_band.dat'.format(prefix) if prefix else 'band.dat'
        if bs.is_metal():
            zero = vr.efermi
        else:
            zero = bs.get_vbm()['energy']

        with open(filename, 'w') as f:
            header = '#k-distance eigenvalue[eV]\n'
            f.write(header)
            for band in bs.bands[Spin.up]:
                for d, e in zip(bs.distance, band):
                    f.write('{:.8f} {:.8f}\n'.format(d, e - zero))
                f.write('\n')

            if bs.is_spin_polarized:
                for band in bs.bands[Spin.down]:
                    for d, e in zip(bs.distance, band):
                        f.write('{:.8f} {:.8f}\n'.format(d, e - zero))
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

    parser.add_argument('-f', '--filenames', default=None, nargs='+',
                        help="one or more vasprun.xml files to plot")
    parser.add_argument('-p', '--prefix', help='prefix for the files generated')
    parser.add_argument('-d', '--directory', help='output directory for files')
    parser.add_argument('-b' '--band-edges', dest='band_edges',
                        action='store_true',
                        help='Highlight the band edges with markers')
    parser.add_argument('--project', default=None, #type=el_orb,
                        help="""Project DOS onto band structure as red, green,
                        and blue contributions. Can project a maximum of 3
                        orbital/elemental contributions. These should be listed
                        using the symbols from the POSCAR, seperated via commas.
                        Specific orbitals can be chosen by adding the orbitals
                        after the element by using a period as the seperator.
                        For example, to project the zinc s and p orbitals, and
                        all the oxygen orbitals, the command would be "--project
                        Zn.s.p,O".""")
    parser.add_argument('--dos', default=None,
                        help="""Path to density of states vasprun.xml.
                        Specifying this option will generate combined DOS/band
                        structure diagrams. The DOS options are more or less
                        the same as for the dosplot command.""")
    parser.add_argument('--elements', type=el_orb, help="""Choose the
                        elements to plot in the DOS. These should be listed
                        using the symbols from the POSCAR and seperated via
                        commas. Specific orbitals can be chosen by adding the
                        orbitals after the element by using a period as the
                        seperator. For example, to plot the carbon s and p, and
                        all the oxygen orbitals, the command would be
                        "--elements C.s.p,O". Must be combined with the --dos
                        option.""")
    parser.add_argument('--orbitals', type=el_orb, help="""Choose the
                        orbitals to split in the DOS. This should be listed as
                        the element (using the symbol from the POSCAR) and the
                        orbitals seperated by a period. For example to plot the
                        oxygen split d orbitals, the command would be
                        "--orbitals O.d". More than one split orbital and
                        element can be added using the notation described for
                        adding more elements. Must be combined with the --dos
                        option.""")
    parser.add_argument('--atoms', type=atoms, help="""Choose which atoms
                        to calculate the DOS for. This should be listed as the
                        element (using the symbol from the POSCAR) and the atoms
                        seperated by a period. For example to plot the oxygen 1,
                        2 and 3 atoms, the command would be "--atoms O.1.2.3".
                        The atom indicies start at 1 (as in the VASP output).
                        You can specify a range to avoid typing all the numbers
                        out, e.g. the previous command can be written "--atoms
                        O.1-3". To select all the atoms of an element just
                        include the element symbol with no numbers after it,
                         e.g. "--atoms Ru" will include all the Ru atoms. If
                        an element is not specified then it will not be
                        included in the DOS. More than one element can be added
                        using the notation described above for adding more
                        elements. Must be combined with the --dos option.""")
    parser.add_argument('--total-only', action='store_true', dest='total_only',
                        help="""Only plot the total DOS. Must be combined with
                        the --dos option""")
    parser.add_argument('--no-total', action='store_false', dest='total',
                        help='Don\'t plot the total DOS')
    parser.add_argument('--legend-cutoff', type=float, default=3,
                        dest='legend_cutoff',
                        help="""Cut-off in %% of total DOS in visible range that
                        determines if a line is given a label. Set to 0 to label
                        all lines. Default is 3 %%. Must be combined with the
                        --dos option.""")
    parser.add_argument('-g', '--gaussian', type=float,
                        help="""Amount of gaussian broadening to apply. Must be
                        combined with the --dos option.""")
    parser.add_argument('--xscale', type=float, default=1,
                        help='Scaling factor for the DOS x axis')
    parser.add_argument('--height', type=float, default=6.0,
                        help='The height of the graph')
    parser.add_argument('--width', type=float, default=6.0,
                        help='The width of the graph')
    parser.add_argument('--ymin', type=float, default=-6.0,
                        help='The minimum energy on the x axis')
    parser.add_argument('--ymax', type=float, default=6.0,
                        help='The maximum energy on the x axis')
    parser.add_argument('--config', type=str, default=None,
                        help='Colour configuration file')
    parser.add_argument('--format', type=str, default='pdf',
                        dest='image_format',
                        help='select image format from pdf, svg, jpg, & png')
    parser.add_argument('--dpi', type=int, default=400,
                        help='pixel density for generated images')
    parser.add_argument('--font', default=None, help='Font to use.')

    args = parser.parse_args()
    logging.basicConfig(filename='vaspy-bandplot.log', level=logging.DEBUG,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    if args.config is None:
        config_path = resource_filename(Requirement.parse('vaspy'),
                                        'conf/orbital_colours.conf')
    else:
        config_path = args.config
    colours = configparser.ConfigParser()
    colours.read(os.path.abspath(config_path))

    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")
    warnings.filterwarnings("ignore", category=UnicodeWarning,
                            module="matplotlib")

    bandplot(filenames=args.filenames, prefix=args.prefix,
             directory=args.directory, vbm_cbm_marker=args.band_edges,
             project=args.project, dos_file=args.dos, elements=args.elements,
             lm_orbitals=args.orbitals, atoms=args.atoms,
             total_only=args.total_only, plot_total=args.total,
             legend_cutoff=args.legend_cutoff, gaussian=args.gaussian,
             height=args.height, width=args.width, ymin=args.ymin,
             ymax=args.ymax, colours=colours, image_format=args.image_format,
             dpi=args.dpi, fonts=[args.font])


if __name__ == "__main__":
    main()
