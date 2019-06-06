# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
A script to plot electronic band structure diagrams.

TODO:
 - Replace the elements and project formats with the dream syntax
"""

from __future__ import unicode_literals
from pkg_resources import Requirement, resource_filename

import os
import sys
import glob
import logging
import argparse
import warnings

from pymatgen.io.vasp.outputs import BSVasprun
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import \
    get_reconstructed_band_structure

import matplotlib as mpl
mpl.use('Agg')

from sumo.io.questaal import QuestaalInit, QuestaalSite, labels_from_syml
from sumo.io.questaal import band_structure as questaal_band_structure
from sumo.plotting.bs_plotter import SBSPlotter
from sumo.plotting.dos_plotter import SDOSPlotter
from sumo.electronic_structure.dos import load_dos
from sumo.cli.dosplot import _atoms, _el_orb

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "July 18, 2017"


def bandplot(filenames=None, code='vasp', prefix=None, directory=None,
             vbm_cbm_marker=False, projection_selection=None, mode='rgb',
             interpolate_factor=4, circle_size=150, dos_file=None,
             cart_coords=False,
             ylabel='Energy (eV)', dos_label=None,
             elements=None, lm_orbitals=None, atoms=None,
             total_only=False, plot_total=True, legend_cutoff=3, gaussian=None,
             height=None, width=None, ymin=-6., ymax=6., colours=None,
             yscale=1, style=None, no_base_style=False,
             image_format='pdf', dpi=400, plt=None, fonts=None):
    """Plot electronic band structure diagrams from vasprun.xml files.

    Args:
        filenames (:obj:`str` or :obj:`list`, optional): Path to input files.
            Vasp:
                Use vasprun.xml or vasprun.xml.gz file.
            Questaal:
                Path to a bnds.ext file. The extension will also be used to
                find site.ext and syml.ext files in the same directory.

            If no filenames are provided, sumo
            will search for vasprun.xml or vasprun.xml.gz files in folders
            named 'split-0*'. Failing that, the code will look for a vasprun in
            the current directory. If a :obj:`list` of vasprun files is
            provided, these will be combined into a single band structure.

        code (:obj:`str`, optional): Calculation type. Default is 'vasp';
            'questaal' also supported (with a reduced feature-set).
        prefix (:obj:`str`, optional): Prefix for file names.
        directory (:obj:`str`, optional): The directory in which to save files.
        vbm_cbm_marker (:obj:`bool`, optional): Plot markers to indicate the
            VBM and CBM locations.
        projection_selection (list): A list of :obj:`tuple` or :obj:`string`
            identifying which elements and orbitals to project on to the
            band structure. These can be specified by both element and
            orbital, for example, the following will project the Bi s, p
            and S p orbitals::

                [('Bi', 's'), ('Bi', 'p'), ('S', 'p')]

            If just the element is specified then all the orbitals of
            that element are combined. For example, to sum all the S
            orbitals::

                [('Bi', 's'), ('Bi', 'p'), 'S']

            You can also choose to sum particular orbitals by supplying a
            :obj:`tuple` of orbitals. For example, to sum the S s, p, and
            d orbitals into a single projection::

                [('Bi', 's'), ('Bi', 'p'), ('S', ('s', 'p', 'd'))]

            If ``mode = 'rgb'``, a maximum of 3 orbital/element
            combinations can be plotted simultaneously (one for red, green
            and blue), otherwise an unlimited number of elements/orbitals
            can be selected.
        mode (:obj:`str`, optional): Type of projected band structure to
            plot. Options are:

                "rgb"
                    The band structure line color depends on the character
                    of the band. Each element/orbital contributes either
                    red, green or blue with the corresponding line colour a
                    mixture of all three colours. This mode only supports
                    up to 3 elements/orbitals combinations. The order of
                    the ``selection`` :obj:`tuple` determines which colour
                    is used for each selection.
                "stacked"
                    The element/orbital contributions are drawn as a
                    series of stacked circles, with the colour depending on
                    the composition of the band. The size of the circles
                    can be scaled using the ``circle_size`` option.
        circle_size (:obj:`float`, optional): The area of the circles used
            when ``mode = 'stacked'``.
        cart_coords (:obj:`bool`, optional): Whether the k-points are read as
            cartesian or reciprocal coordinates. This is only required for
            Questaal output; Vasp output is less ambiguous. Defaults to
            ``False`` (fractional coordinates).
        dos_file (:obj:'str', optional): Path to vasprun.xml file from which to
            read the density of states information. If set, the density of
            states will be plotted alongside the bandstructure.
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
        total_only (:obj:`bool`, optional): Only extract the total density of
            states. Defaults to ``False``.
        plot_total (:obj:`bool`, optional): Plot the total density of states.
            Defaults to ``True``.
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
        ymin (:obj:`float`, optional): The minimum energy on the y-axis.
        ymax (:obj:`float`, optional): The maximum energy on the y-axis.
        style (:obj:`list` or :obj:`str`, optional): (List of) matplotlib style
            specifications, to be composed on top of Sumo base style.
        no_base_style (:obj:`bool`, optional): Prevent use of sumo base style.
            This can make alternative styles behave more predictably.
        colours (:obj:`dict`, optional): Use custom colours for specific
            element and orbital combinations. Specified as a :obj:`dict` of
            :obj:`dict` of the colours. For example::

                {
                    'Sn': {'s': 'r', 'p': 'b'},
                    'O': {'s': '#000000'}
                }

            The colour can be a hex code, series of rgb value, or any other
            format supported by matplotlib.
        yscale (:obj:`float`, optional): Scaling factor for the y-axis.
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
        If ``plt`` set then the ``plt`` object will be returned. Otherwise, the
        method will return a :obj:`list` of filenames written to disk.
    """
    if not filenames:
        filenames = find_vasprun_files()
    elif isinstance(filenames, str):
        filenames = [filenames]

    # only load the orbital projects if we definitely need them
    parse_projected = True if projection_selection else False

    # now load all the band structure data and combine using the
    # get_reconstructed_band_structure function from pymatgen
    bandstructures = []
    if code == 'vasp':
        for vr_file in filenames:
            vr = BSVasprun(vr_file, parse_projected_eigen=parse_projected)
            bs = vr.get_band_structure(line_mode=True)
            bandstructures.append(bs)
        bs = get_reconstructed_band_structure(bandstructures)
    elif code == 'questaal':
        bnds_file = filenames[0]
        ext = bnds_file.split('.')[-1]
        bnds_folder = os.path.join(bnds_file, os.path.pardir)

        site_file = os.path.abspath(
            os.path.join(bnds_folder, 'site.{}'.format(ext)))

        if os.path.isfile(site_file):
            logging.info('site file found, reading lattice...')
            site_data = QuestaalSite.from_file(site_file)
            #bnds_lattice = QuestaalInit.from_file(init_file).structure.lattice
            bnds_lattice = site_data.structure.lattice
            alat = site_data.alat
        else:
            raise IOError('Site file {} not found: '
                          'needed to determine lattice'.format(site_file))

        syml_file = os.path.abspath(
            os.path.join(bnds_folder, 'syml.{}'.format(ext)))
        if os.path.isfile(syml_file):
            logging.info('syml file found, reading special-point labels...')
            bnds_labels = labels_from_syml(syml_file)
        else:
            logging.info('syml file not found, band structure lacks labels')
            bnds_labels = {}

        bs = questaal_band_structure(bnds_file, bnds_lattice, alat=alat,
                                     labels=bnds_labels,
                                     coords_are_cartesian=cart_coords)

    # currently not supported as it is a pain to make subplots within subplots,
    # although need to check this is still the case
    if 'split' in mode and dos_file:
        logging.error('ERROR: Plotting split projected band structure with DOS'
                      ' not supported.\nPlease use --projected-rgb or '
                      '--projected-stacked options.')
        sys.exit()

    if (projection_selection and mode == 'rgb'
            and len(projection_selection) > 3):
        logging.error('ERROR: RGB projected band structure only '
                      'supports up to 3 elements/orbitals.'
                      '\nUse alternative --mode setting.')
        sys.exit()

    # don't save if pyplot object provided
    save_files = False if plt else True

    dos_plotter = None
    dos_opts = None
    if dos_file:
        dos, pdos = load_dos(dos_file, elements, lm_orbitals, atoms, gaussian,
                             total_only)
        dos_plotter = SDOSPlotter(dos, pdos)
        dos_opts = {'plot_total': plot_total, 'legend_cutoff': legend_cutoff,
                    'colours': colours, 'yscale': yscale}

    plotter = SBSPlotter(bs)
    if projection_selection:
        plt = plotter.get_projected_plot(
            projection_selection, mode=mode,
            interpolate_factor=interpolate_factor, circle_size=circle_size,
            zero_to_efermi=True, ymin=ymin, ymax=ymax, height=height,
            width=width, vbm_cbm_marker=vbm_cbm_marker, ylabel=ylabel,
            plt=plt, dos_plotter=dos_plotter, dos_options=dos_opts,
            dos_label=dos_label, fonts=fonts, style=style,
            no_base_style=no_base_style)
    else:
        plt = plotter.get_plot(
            zero_to_efermi=True, ymin=ymin, ymax=ymax, height=height,
            width=width, vbm_cbm_marker=vbm_cbm_marker, ylabel=ylabel,
            plt=plt, dos_plotter=dos_plotter, dos_options=dos_opts,
            dos_label=dos_label, fonts=fonts, style=style,
            no_base_style=no_base_style)

    if save_files:
        basename = 'band.{}'.format(image_format)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename
        if directory:
            filename = os.path.join(directory, filename)
        plt.savefig(filename, format=image_format, dpi=dpi,
                    bbox_inches='tight')

        written = [filename]
        written += save_data_files(bs, prefix=prefix,
                                   directory=directory)
        return written

    else:
        return plt


def find_vasprun_files():
    """Search for vasprun files from the current directory.

    The precedence order for file locations is:

      1. First search for folders named: 'split-0*'
      2. Else, look in the current directory.

    The split folder names should always be zero based, therefore easily
    sortable.
    """
    folders = glob.glob('split-*')
    folders = sorted(folders) if folders else ['.']

    filenames = []
    for fol in folders:
        vr_file = os.path.join(fol, 'vasprun.xml')
        vr_file_gz = os.path.join(fol, 'vasprun.xml.gz')

        if os.path.exists(vr_file):
            filenames.append(vr_file)
        elif os.path.exists(vr_file_gz):
            filenames.append(vr_file_gz)
        else:
            logging.error('ERROR: No vasprun.xml found in {}!'.format(fol))
            sys.exit()

    return filenames


def save_data_files(bs, prefix=None, directory=None):
    """Write the band structure data files to disk.

    Args:
        bs (`BandStructureSymmLine`): Calculated band structure.
        prefix (`str`, optional): Prefix for data file.
        directory (`str`, optional): Directory in which to save the data.

    Returns:
        The filename of the written data file.
    """
    filename = '{}_band.dat'.format(prefix) if prefix else 'band.dat'
    directory = directory if directory else '.'
    filename = os.path.join(directory, filename)

    if bs.is_metal():
        zero = bs.efermi
    else:
        zero = bs.get_vbm()['energy']

    with open(filename, 'w') as f:
        header = '#k-distance eigenvalue[eV]\n'
        f.write(header)

        # write the spin up eigenvalues
        for band in bs.bands[Spin.up]:
            for d, e in zip(bs.distance, band):
                f.write('{:.8f} {:.8f}\n'.format(d, e - zero))
            f.write('\n')

        # calculation is spin polarised, write spin down bands at end of file
        if bs.is_spin_polarized:
            for band in bs.bands[Spin.down]:
                for d, e in zip(bs.distance, band):
                    f.write('{:.8f} {:.8f}\n'.format(d, e - zero))
                f.write('\n')
    return filename


def _el_orb_tuple(string):
    """Parse the element and orbital argument strings.

    The presence of an element without any orbitals means that we want to plot
    all of its orbitals.

    Args:
        string (`str`): The selected elements and orbitals in in the form:
            `"Sn.s.p,O"`.

    Returns:
        A list of tuples specifying which elements/orbitals to plot. The output
        for the above example would be:

            `[('Sn', ('s', 'p')), 'O']`
    """
    el_orbs = []
    for split in string.split(','):
        splits = split.split('.')
        el = splits[0]
        if len(splits) == 1:
            el_orbs.append(el)
        else:
            el_orbs.append((el, tuple(splits[1:])))
    return el_orbs


def _get_parser():
    parser = argparse.ArgumentParser(description="""
    bandplot is a script to produce publication-ready band
    structure diagrams""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filenames', default=None, nargs='+',
                        metavar='F',
                        help="one or more vasprun.xml files to plot")
    parser.add_argument('-c', '--code', default='vasp',
                        help='Electronic structure code (default: vasp).'
                             '"questaal" also supported.')
    parser.add_argument('-p', '--prefix', metavar='P',
                        help='prefix for the files generated')
    parser.add_argument('-d', '--directory', metavar='D',
                        help='output directory for files')
    parser.add_argument('-b', '--band-edges', dest='band_edges',
                        action='store_true',
                        help='highlight the band edges with markers')
    parser.add_argument('--project', default=None, metavar='S',
                        type=_el_orb_tuple, dest='projection_selection',
                        help=('select which orbitals to project onto the band '
                              'structure (e.g. "Zn.s,Zn.p,O")'))
    parser.add_argument('--mode', default='rgb', type=str,
                        help=('mode for orbital projections (options: rgb, '
                              'stacked)'))
    parser.add_argument('--interpolate-factor', type=int, default=4,
                        dest='interpolate_factor', metavar='N',
                        help=('interpolate factor for band structure '
                              'projections (default: 4)'))
    parser.add_argument('--cartesian', action='store_true',
                        help='Read cartesian k-point coordinates. This is only'
                             ' necessary for some Questaal calculations; Vasp '
                             'outputs are less ambiguous and this option will '
                             'be ignored if --code=vasp.')
    parser.add_argument('--circle-size', type=int, default=150,
                        dest='circle_size', metavar='S',
                        help=('circle size for "stacked" projections '
                              '(default: 150)'))
    parser.add_argument('--ylabel', type=str, default='Energy (eV)',
                        help='y-axis (i.e. energy) label/units')
    parser.add_argument('--dos-label', type=str, dest='dos_label',
                        default=None,
                        help='Axis label for DOS if included')
    parser.add_argument('--dos', default=None,
                        help='path to density of states vasprun.xml')
    parser.add_argument('--elements', type=_el_orb, metavar='E',
                        help='elemental orbitals to plot (e.g. "C.s.p,O")')
    parser.add_argument('--orbitals', type=_el_orb, metavar='O',
                        help=('orbitals to split into lm-decomposed '
                              'contributions (e.g. "Ru.d")'))
    parser.add_argument('--atoms', type=_atoms, metavar='A',
                        help=('atoms to include (e.g. "O.1.2.3,Ru.1.2.3")'))
    parser.add_argument('--total-only', action='store_true', dest='total_only',
                        help='only plot the total density of states')
    parser.add_argument('--no-total', action='store_false', dest='total',
                        help='don\'t plot the total density of states')
    parser.add_argument('--legend-cutoff', type=float, default=3,
                        dest='legend_cutoff', metavar='C',
                        help=('cut-off in %% of total DOS that determines if'
                              ' a line is given a label (default: 3)'))
    parser.add_argument('-g', '--gaussian', type=float, metavar='G',
                        help='standard deviation of DOS gaussian broadening')
    parser.add_argument('--scale', type=float, default=1,
                        help='scaling factor for the density of states')
    parser.add_argument('--height', type=float, default=None,
                        help='height of the graph')
    parser.add_argument('--width', type=float, default=None,
                        help='width of the graph')
    parser.add_argument('--ymin', type=float, default=-6.,
                        help='minimum energy on the y-axis')
    parser.add_argument('--ymax', type=float, default=6.,
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
    parser.add_argument('--dpi', type=int, default=400,
                        help='pixel density for image file')
    parser.add_argument('--font', default=None, help='font to use')
    return parser


def main():
    args = _get_parser().parse_args()
    logging.basicConfig(filename='sumo-bandplot.log', level=logging.INFO,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    if args.config is None:
        config_path = resource_filename(Requirement.parse('sumo'),
                                        'sumo.plotting/orbital_colours.conf')
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

    bandplot(filenames=args.filenames, code=args.code, prefix=args.prefix,
             directory=args.directory, vbm_cbm_marker=args.band_edges,
             projection_selection=args.projection_selection, mode=args.mode,
             interpolate_factor=args.interpolate_factor,
             cart_coords=args.cartesian,
             circle_size=args.circle_size, yscale=args.scale,
             ylabel=args.ylabel, dos_label=args.dos_label,
             dos_file=args.dos, elements=args.elements,
             lm_orbitals=args.orbitals, atoms=args.atoms,
             total_only=args.total_only, plot_total=args.total,
             legend_cutoff=args.legend_cutoff, gaussian=args.gaussian,
             height=args.height, width=args.width, ymin=args.ymin,
             style=args.style, no_base_style=args.no_base_style,
             ymax=args.ymax, colours=colours, image_format=args.image_format,
             dpi=args.dpi, fonts=args.font)


if __name__ == "__main__":
    main()
