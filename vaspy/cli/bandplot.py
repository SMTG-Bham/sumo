# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals
from pkg_resources import Requirement, resource_filename

import os
import sys
import glob
import logging
import argparse
import warnings

import matplotlib as mpl
mpl.use('Agg')

from vaspy.plotting.bs_plotter import VBSPlotter
from vaspy.plotting.dos_plotter import VDOSPlotter
from vaspy.electronic_structure.dos import load_dos
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

# TODO:
#  - Replace the elements and project formats with the dream syntax


def bandplot(filenames=None, prefix=None, directory=None, vbm_cbm_marker=False,
             project_split=None, project_rgb=None, project_stacked=None,
             interpolate_factor=4, circle_size=150, dos_file=None,
             elements=None, lm_orbitals=None, atoms=None,
             total_only=False, plot_total=True, legend_cutoff=3, gaussian=None,
             height=6., width=6., ymin=-6., ymax=6., colours=None, yscale=1,
             image_format='pdf', dpi=400, plt=None, fonts=None):
    """Plot electronic band structure diagrams from vasprun.xml files.

    Args:
        filenames (`str`, `list`, optional): Names of vasprun.xml files to use
            in the band structure (files can be gziped). If no filenames are
            provided, the code will search for vasprun.xml or vasprun.xml.gz
            files in folders named 'split-0*'. Failing that, the code will look
            for a vasprun in the current directory. If a `list` of vasprun
            files is provided, these will be combined into a single band
            structure.
        prefix (`str`, optional): Prefix for files generated.
        directory (`str`, optional): Directory in which to save files.
        vbm_cbm_marker (`bool`, optional): Mark the valence band maxima and
            conduction band minima using coloured circles.
        project_split: WIP
        project_rgb (`list`, optional): Which orbitals to project onto the band
            structure using red, green, and blue. Should be specified as a
            `list` of `tuple` or `str`, with the order of the `tuple`s defining
            the colour used. For example, to plot the Sn s, p, and d orbitals
            as red, green, and blue, respectively, `project_rgb` should be:

                `[('Sn', 's'), ('Sn', 'p'), ('Sn', 'd')]`

            Multiple orbitals can be summed together into a single projection.
            For example, the following will  plot the Sn s and p orbitals as
            red, Sn d orbitals as green, and oxygen p orbitals as blue:

                `[('Sn', ('s', 'p')), ('Sn', 'd'), ('O', 'p')]`

            Just suppling the element name is a shorthanden for summing all of
            its orbital contributions. For example, the following will plot the
            sum of all Sn orbitals as red, all O orbitals as blue, and all Cs
            orbitals as green:

                `['Sn', 'O', 'Cs']`

            A maximum of 3 orbitals can be projected simultaneously.
        project_stacked (`list`, optional): Write this when projection cli is
            finished.
        interpolate_factor (`int`, optional): The factor by which to
            interpolate the band structure (neccessary to make smooth lines for
            projected plots). A larger number indicates greater interpolation.
        circle_size (`float`, optional): The size of the circles used in the
            project_stacked plotting mode.
        dos_file ('str', optional): vasprun.xml file from which to read the
            density of states information. If set, the density of states will
            be plotted alongside the bandstructure.
        elements (`dict`, optional): The elements and orbitals to plot in the
            density of states. Should be provided as a `dict` with the keys as
            the element names and corresponding values as a `tuple` of orbitals
            to plot. For example, the following would plot the Bi s, px, py and
            d orbitals:

                `{'Bi': ('s', 'px', 'py', 'd')}`.

            If an element is included with an empty `tuple`, all orbitals for
            that species will be plotted. If `elements` is not set or set to
            `None`, all elements for all species will be considered.
        lm_orbitals (`dict`, optional): The orbitals to decompose into their lm
            contributions (e.g. p -> px, py, pz). Should be provided as a
            `dict`, with the elements names as keys and a `tuple` of orbitals
            as the corresponding values. For example, the following would be
            used to decompose the oxygen p and d orbitals:

                `{'O': ('p', 'd')}'

        atoms (`dict`, optional): Which atomic sites to plot the density of
            states for. Should be provided as a `dict`, with the element names
            as keys and a `tuple` of `int` specifiying the atomic indicies as
            the corresponding values. The elemental projected density of states
            will be summed only over the atom inidices specified. If an element
            is included with an empty `tuple`, then all sites for that element
            will be included. The indices are 0 based for each element
            specified in the POSCAR. For example, the following will calculate
            the denisty of states for the first 4 Sn atoms and all O atoms in
            the structure:

                `{'Sn': (1, 2, 3, 4), 'O': (, )}`

            If `atoms` is not set or set to `None` then all atomic sites for
            all elements will be considered.
        total_only (`bool`, optional): Only plot the total density of states.
        plot_total (:obj:`bool`, optional): Plot the total density of states.
            Defaults to ``True``.
        legend_cutoff (:obj:`float`, optional): The cut-off (in % of the
            maximum density of states within the plotting range) for an
            elemental orbital to be labelled in the legend. This prevents the
            legend from containing labels for orbitals that have very little
            contribution in the plotting range.
        gaussian (`float`, optional): Broaden the density of states using
            convolution with a gaussian function. This parameter controls the
            sigma or smearing width of the gaussian.
        height (`float`, optional): The height of the plot in inches.
        width (`float`, optional): The width of the plot in inches.
        xmin (`float`, optional): The minimum energy to plot.
        xmax (`float`, optional): The maximum energy to plot.
        colours (`dict`, optional): Colours to use when plotting elemental
            density of states. Should be provided as a `dict`, where the key is
            the element name and the corresponding value is a `dict` of
            orbitals and their colour. The colour can be any matplotlib
            supported colour identifier, e.g. hex, rgb, or name. For example,
            the following will set the O p orbitals to red and the Sn s
            orbitals to green.

                `{'Sn': {'s': 'r'}, 'O': {'p': 'g'}}`

            If an orbital colour is not specified, the code will select a
            colour from a list of 21 visually distinct colours.
        yscale (`float`, optional): Scaling factor for the y-axis.
        image_format (`str`, optional): The image file format. Can be any
            format supported by matplot, including: png, jpg, pdf, and svg.
        dpi (`int`, optional): The dots-per-inch (pixel density) for the image.
        plt (`matplotlib.pyplot`, optional): Matplotlib object to use for
            plotting. If plt is set then no files will be written.
        fonts (`list`, optional): A `list` of fonts to try and use. Preference
            will be given to the fonts at he beginning of the list.

    Returns:
        If `plt` set then the `plt object will be returned. Otherwise, the
        method will return a `list` of filenames written to disk.
    """
    if not filenames:
        filenames = find_vasprun_files()
    elif type(filenames) == str:
        filenames = [filenames]

    # only laod the orbital proejcts if we definitely need them
    parse_projected = True if (project_split or project_rgb
                               or project_stacked) else False

    # now load all the vaspruns and combine them together using the
    # get_reconstructed_band_structure function from pymatgen
    bandstructures = []
    for vr_file in filenames:
        vr = BSVasprun(vr_file, parse_projected_eigen=parse_projected)
        bs = vr.get_band_structure(line_mode=True)
        bandstructures.append(bs)
    bs = get_reconstructed_band_structure(bandstructures)

    # currently not supported as it is a pain to make subplots within subplots,
    # although need to check this is still the case
    if project_split and dos_file:
        logging.error('ERROR: Plotting split projected band structure with DOS'
                      ' not supported.\nPlease use --projected-rgb or '
                      '--projected-stacked options.')
        sys.exit()

    if project_rgb and len(project_rgb) > 3:
        logging.error('ERROR: Plotting RGB projected band structure only '
                      'supports up to 3 elements/orbitals.')
        sys.exit()

    # don't save if pyplot object provided
    save_files = False if plt else True

    dos_plotter = None
    dos_opts = None
    if dos_file:
        dos, pdos = load_dos(dos_file, elements, lm_orbitals, atoms, gaussian,
                             total_only)
        dos_plotter = VDOSPlotter(dos, pdos)
        dos_opts = {'plot_total': plot_total, 'legend_cutoff': legend_cutoff,
                    'colours': colours, 'yscale': yscale}

    plotter = VBSPlotter(bs)
    if project_rgb or project_split or project_stacked:

        # projected plotter logic could probably be improved
        if project_rgb:
            mode = 'rgb'
            selection = project_rgb
        elif project_split:
            mode = 'split'
            selection = project_split
            raise NotImplementedError('projected band structure plotting not '
                                      'yet supported')
        elif project_stacked:
            selection = project_stacked
            mode = 'stacked'

        plt = plotter.get_projected_plot(selection, mode=mode,
                                         interpolate_factor=interpolate_factor,
                                         circle_size=circle_size,
                                         zero_to_efermi=True, ymin=ymin,
                                         ymax=ymax, height=height, width=width,
                                         vbm_cbm_marker=vbm_cbm_marker,
                                         plt=plt, dos_plotter=dos_plotter,
                                         dos_options=dos_opts, fonts=fonts)
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
        plt.savefig(filename, format=image_format, dpi=dpi,
                    bbox_inches='tight')

        written = [filename]
        written += save_data_files(bs, prefix=prefix, directory=directory)
        return written

    else:
        return plt


def find_vasprun_files():
    """Search for vasprun files from the current directory.

    The precedence order for file locations is:

      1. First search for folders named: 'split-0*'
      2. Else, look in the current directory.

    The split folder names should always be zerop based, therefore easily
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


def save_data_files(vr, bs, prefix=None, directory=None):
    """Write the band structure data files to disk.

    Args:
        vs (`Vasprun`): Pymatgen `Vasprun` object.
        bs (`BandStructureSymmLine`): Calculated band structure.
        prefix (`str`, optional): Prefix for data file.
        directory (`str`, optional): Directory in which to save the data.

    Returns:
        The filename of the written data file.
    """
    filename = '{}_band.dat'.format(prefix) if prefix else 'band.dat'
    filename = os.path.join(directory, filename)

    if bs.is_metal():
        zero = vr.efermi
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


def el_orb_tuple(string):
    """Parse the element and orbital argument strings.

    The presence of an element without any orbitals means that we want to plot
    all of its orbitals.

    Args:
        string (`str`): The selected eleemtns and orbitals in in the form:
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


def main():
    parser = argparse.ArgumentParser(description="""
    bandplot is a convenient script to help make publication ready band
    structure diagrams.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filenames', default=None, nargs='+',
                        help="one or more vasprun.xml files to plot")
    parser.add_argument('-p', '--prefix', help='prefix for generated files')
    parser.add_argument('-d', '--directory', help='output directory for files')
    parser.add_argument('-b' '--band-edges', dest='band_edges',
                        action='store_true',
                        help='Highlight the band edges with markers')
    parser.add_argument('--project-rgb', default=None, type=el_orb_tuple,
                        dest='project_rgb',
                        help="""Project orbital contributions onto band
                        structure as red, green, and blue lines. Can project
                        a maximum of 3
                        orbital/elemental contributions. These should be listed
                        using the symbols from the POSCAR, seperated via
                        commas. Specific orbitals can be chosen by adding the
                        orbitals after the element by using a period as the
                        seperator. For example, to project the zinc s as red,
                        zinc p as green, and sum all oxygen atoms as blue,
                        the command would be "--project-rgb Zn.s,Zn.p,O".""")
    parser.add_argument('--project-stacked', default=None, type=el_orb_tuple,
                        dest='project_stacked',
                        help="""Project orbtal contributions onto band
                        structure as a series of coloured circles. These should
                        be listed using the symbols from the POSCAR, seperated
                        via commas. Specific orbitals can be chosen by adding
                        the orbitals after the element by using a period as the
                        seperator. For example, to project the zinc s as red,
                        zinc p as green, and sum all oxygen atoms as blue,
                        the command would be "--project-rgb Zn.s,Zn.p,O".""")
    parser.add_argument('-i', '--interpolate-factor', type=int, default=4,
                        dest='interpolate_factor',
                        help="""For projected band structures only: Interpolate
                        the band structure and projections. A larger factor
                        indicates increased interpolate. Default is 4.""")
    parser.add_argument('--circle-size', type=int, default=150,
                        dest='circle_size',
                        help="""For stacked projected band structures only:
                        Maximum size of circles drawn. Default is 150.""")
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
                        element (using the symbol from the POSCAR) and the
                        atoms seperated by a period. For example to plot the
                        oxygen 1, 2 and 3 atoms, the command would be "--atoms
                        O.1.2.3". The atom indicies start at 1 (as in the VASP
                        output). You can specify a range to avoid typing all
                        the numbers out, e.g. the previous command can be
                        written "--atoms O.1-3".
                        To select all the atoms of an element just
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
                        help="""Cut-off in %% of total DOS in plotting range
                        that determines if a line is given a label. Set to 0 to
                        label all lines. Default is 3 %%. Must be combined with
                        the --dos option.""")
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
             project_rgb=args.project_rgb,
             project_stacked=args.project_stacked,
             interpolate_factor=args.interpolate_factor,
             circle_size=args.circle_size,
             dos_file=args.dos, elements=args.elements,
             lm_orbitals=args.orbitals, atoms=args.atoms,
             total_only=args.total_only, plot_total=args.total,
             legend_cutoff=args.legend_cutoff, gaussian=args.gaussian,
             height=args.height, width=args.width, ymin=args.ymin,
             ymax=args.ymax, colours=colours, image_format=args.image_format,
             dpi=args.dpi, fonts=[args.font])


if __name__ == "__main__":
    main()
