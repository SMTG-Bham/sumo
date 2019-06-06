# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
A script to calculate and plot optical spectra from ab initio calculations.
"""

import os
from glob import glob
import sys
import logging
import warnings
import argparse
from collections import OrderedDict

import matplotlib as mpl
mpl.use('Agg')

from pymatgen.io.vasp import Vasprun
from pymatgen.util.string import latexify
from sumo.io import questaal

from sumo.plotting.optics_plotter import SOpticsPlotter
from sumo.electronic_structure.optics import (broaden_eps,
                                              calculate_dielectric_properties,
                                              write_files)

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Jan 10, 2018"


def optplot(modes=('absorption',), filenames=None, codes='vasp',
            prefix=None, directory=None,
            gaussian=None, band_gaps=None, labels=None, average=True, height=6,
            width=6, xmin=0, xmax=None, ymin=0, ymax=1e5, colours=None,
            style=None, no_base_style=None,
            image_format='pdf', dpi=400, plt=None, fonts=None):
    """A script to plot optical absorption spectra from VASP calculations.

    Args:
        modes (:obj:`list` or :obj:`tuple`):
            Ordered list of :obj:`str` determining properties to plot.
            Accepted options are 'absorption' (default), 'eps', 'eps-real',
                'eps-im', 'n', 'n-real', 'n-im', 'loss' (equivalent to n-im).
        filenames (:obj:`str` or :obj:`list`, optional): Path to data file.
            For VASP this is a *vasprun.xml* file (can be gzipped); for
            Questaal the *opt.ext* file from *lmf* or *eps_BSE.out* from
            *bethesalpeter* may be used.
            Alternatively, a list of paths can be
            provided, in which case the absorption spectra for each will be
            plotted concurrently.
        codes (:obj:`str` or :obj:`list`, optional): Original
            calculator. Accepted values are 'vasp' and 'questaal'. Items should
            correspond to filenames.
        prefix (:obj:`str`, optional): Prefix for file names.
        directory (:obj:`str`, optional): The directory in which to save files.
        gaussian (:obj:`float`): Standard deviation for gaussian broadening.
        band_gaps (:obj:`float`, :obj:`str` or :obj:`list`, optional): The band
            gap as a :obj:`float`, plotted as a dashed line. If plotting
            multiple spectra then a :obj:`list` of band gaps can be provided.
            Band gaps can be provided as a floating-point number or as a path
            to a *vasprun.xml* file. To skip over a line, set its bandgap to
            zero or a negative number to place it outside the visible range.
        labels (:obj:`str` or :obj:`list`): A label to identify the spectra.
            If plotting multiple spectra then a :obj:`list` of labels can
            be provided.
        average (:obj:`bool`, optional): Average the dielectric response across
            all lattice directions. Defaults to ``True``.
        height (:obj:`float`, optional): The height of the plot.
        width (:obj:`float`, optional): The width of the plot.
        xmin (:obj:`float`, optional): The minimum energy on the x-axis.
        xmax (:obj:`float`, optional): The maximum energy on the x-axis.
        ymin (:obj:`float`, optional): The minimum absorption intensity on the
            y-axis.
        ymax (:obj:`float`, optional): The maximum absorption intensity on the
            y-axis.
        colours (:obj:`list`, optional): A :obj:`list` of colours to use in the
            plot. The colours can be specified as a hex code, set of rgb
            values, or any other format supported by matplotlib.
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

    # Don't write files if this is being done to manipulate existing plt
    save_files = False if plt else True

    ##### BUILD LIST OF FILES AUTOMATICALLY IF NECESSARY #####

    if codes == 'vasp':
        if not filenames:
            if os.path.exists('vasprun.xml'):
                filenames = ['vasprun.xml']
            elif os.path.exists('vasprun.xml.gz'):
                filenames = ['vasprun.xml.gz']
            else:
                logging.error('ERROR: No vasprun.xml found!')
                sys.exit()

    elif codes == 'questaal':
        if not filenames:
            if len(glob('opt.*')) > 0:
                filenames = glob('opt.*')
                if len(filenames) == 1:
                    logging.info("Found optics file: " + filenames[0])
                else:
                    logging.info("Found optics files: " + ", ".join(filenames))

    if isinstance(filenames, str):
        filenames = [filenames]

    if isinstance(codes, str):
        codes = [codes] * len(filenames)
    elif len(codes) == 1:
        codes = list(codes) * len(filenames)

    #### ITERATE OVER FILES READING DIELECTRIC DATA ####

    dielectrics = []
    auto_labels = []
    auto_band_gaps = []
    for i, (filename, code) in enumerate(zip(filenames, codes)):
        if code == 'vasp':
            vr = Vasprun(filename)
            dielectrics.append(vr.dielectric)

            auto_labels.append(
                latexify(vr.final_structure.composition.reduced_formula).
                replace('$_', '$_\mathregular'))

            if isinstance(band_gaps, list) and not band_gaps:
                # band_gaps = [], auto band gap requested
                auto_band_gaps.append(
                    vr.get_band_structure().get_band_gap()['energy'])
            else:
                auto_band_gaps.append(None)

        elif code == 'questaal':
            if not save_files:
                out_filename = None
            elif len(filenames) == 1:
                out_filename = 'dielectric.dat'
            else:
                out_filename = 'dielectric_{0}.dat'.format(i + 1)

            dielectrics.append(
                questaal.dielectric_from_file(filename, out_filename))

            auto_band_gaps.append(None)
            auto_labels.append(filename.split('.')[-1])
            if isinstance(band_gaps, list) and not band_gaps:
                logging.info('Bandgap requested but not supported for Questaal'
                             ' file {}: skipping...'.format(filename))

        else:
            raise Exception('Code selection "{}" not recognised'.format(code))

    if not labels and len(filenames) > 1:
        labels = auto_labels

    #### PROCESS DIELECTRIC DATA: BROADENING AND DERIVED PROPERTIES ####

    if gaussian:
        dielectrics = [broaden_eps(d, gaussian)
                       for d in dielectrics]

    # initialize spectrum data ready to append from each dataset
    abs_data = OrderedDict()

    for mode in modes:
        abs_data.update({mode: []})

    # for each calculation, get all required properties and append to data
    for d in dielectrics:
        for mode, spectrum in calculate_dielectric_properties(
                d, set(modes), average=average).items():
            abs_data[mode].append(spectrum)

    if isinstance(band_gaps, list) and not band_gaps:
        # empty list therefore use bandgaps collected from vasprun files
        band_gaps = auto_band_gaps
    elif isinstance(band_gaps, list):
        # list containing filenames and/or values: mutate the list in-place
        for i, item in enumerate(band_gaps):
            if item is None:
                pass
            elif _floatable(item):
                band_gaps[i] = float(item)
            elif 'vasprun' in item:
                band_gaps[i] = (
                    Vasprun(item).get_band_structure().get_band_gap()['energy']
                    )
            else:
                raise ValueError('Format not recognised for auto bandgap: '
                                 '{}.'.format(item))

    plotter = SOpticsPlotter(abs_data, band_gap=band_gaps, label=labels)
    plt = plotter.get_plot(width=width, height=height, xmin=xmin,
                           xmax=xmax, ymin=ymin, ymax=ymax,
                           colours=colours, dpi=dpi, plt=plt, fonts=fonts,
                           style=style, no_base_style=no_base_style)

    if save_files:
        basename = 'absorption'
        if prefix:
            basename = '{}_{}'.format(prefix, basename)
        image_filename = '{}.{}'.format(basename, image_format)

        if directory:
            image_filename = os.path.join(directory, image_filename)
        plt.savefig(image_filename, format=image_format, dpi=dpi)
        for mode, data in abs_data.items():
            basename = 'absorption' if mode == 'abs' else mode
            write_files(data, basename=basename,
                        prefix=prefix, directory=directory)
    else:
        return plt

def _floatable(item):
    """Check if an item can be intepreted with float()"""
    try:
        float(item)
        return True
    except ValueError:
        return False

def _get_parser():
    parser = argparse.ArgumentParser(description="""
    optplot is a script to produce optical absorption spectra diagrams""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('mode', type=str, nargs='*', default='absorption',
                        metavar='M',
                        choices={'absorption', 'loss', 'eps_real', 'eps_imag',
                                 'n_real', 'n_imag'},
                        help='Optical properties to plot. Multiple choices '
                             ' will be displayed as subplots. Accepted values:'
                             ' "absorption" (optical absorption over distance)'
                             ', "loss" (energy-loss function -Im(1/eps)), '
                             '"eps_real" and "eps_imag" (real and imaginary '
                             'parts of the dielectric function), '
                             '"n_real" (real part of complex refractive index)'
                             '"n_imag" (imaginary part of RI, also known as '
                             'the extinction coefficient kappa.)')
    parser.add_argument('-f', '--filenames', metavar='F',
                        help='path to one or more vasprun.xml files',
                        default=None, nargs='+')
    parser.add_argument('-p', '--prefix', metavar='P',
                        help='prefix for the files generated')
    parser.add_argument('-d', '--directory', metavar='D',
                        help='output directory for files')
    parser.add_argument('-c', '--code', metavar='C', default='vasp', nargs='+',
                        help=('Original calculator. Accepted values are '
                              '"vasp" and "questaal".'))
    parser.add_argument('-g', '--gaussian', type=float, metavar='G',
                        help='standard deviation of gaussian broadening')
    parser.add_argument('-b', '--bandgaps', nargs='*', metavar='E',
                        help=('indicate the fundamental band gap (options: '
                              'nothing, vasprun.xml file, or float). A '
                              'sequence of files and values may be provided, '
                              'corresponding to the optical data files. '
                              'To skip a line, set a value outside the plot '
                              'range (e.g. -1).'))
    parser.add_argument('-l', '--labels', nargs='+', metavar='L',
                        help='labels for the absorption specta')
    parser.add_argument('-a', '--anisotropic', action='store_false',
                        help='separate spectra into to x, y, and z directions')
    parser.add_argument('--height', type=float, default=None,
                        help='height of the graph')
    parser.add_argument('--width', type=float, default=None,
                        help='width of the graph')
    parser.add_argument('--xmin', type=float, default=0.,
                        help='minimum energy on the x-axis')
    parser.add_argument('--xmax', type=float, default=None,
                        help='maximum energy on the x-axis')
    parser.add_argument('--ymin', type=str, default=['auto'], nargs='+',
                        help='minimum intensity on the y-axis; may specify '
                             'multiple values if plotting more than one axis. '
                             'Use "auto" or "_" for automatic value.')
    parser.add_argument('--ymax', type=str, default=['auto'], nargs='+',
                        help='maximum intensity on the y-axis; may specify'
                             'multiple values if plotting more than one axis. '
                             'Use "auto" or "_" for automatic value.')
    parser.add_argument('--style', type=str, nargs='+', default=None,
                        help='matplotlib style specifications')
    parser.add_argument('--no-base-style', action='store_true',
                        dest='no_base_style',
                        help='prevent use of sumo base style')
    parser.add_argument('--format', type=str, default='pdf',
                        dest='image_format', metavar='FORMAT',
                        help='image file format (options: pdf, svg, jpg, png)')
    parser.add_argument('--dpi', type=int, default=400,
                        help='pixel density for image file')
    parser.add_argument('--font', default=None, help='font to use')
    return parser


def main():
    args = _get_parser().parse_args()

    logging.basicConfig(filename='sumo-optplot.log', level=logging.INFO,
                        filemode='w', format='%(message)s')
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger('').addHandler(console)

    warnings.filterwarnings("ignore", category=UserWarning,
                            module="matplotlib")
    warnings.filterwarnings("ignore", category=UnicodeWarning,
                            module="matplotlib")
    warnings.filterwarnings("ignore", category=UserWarning,
                            module="pymatgen")

    # Wrap mode into list if necessary
    if not isinstance(args.mode, list):
        args.mode = [args.mode]

    # Replace text placeholders with preferred Python representation: None
    ymin = [None if (x.lower() in ('auto', '_')) else float(x)
            for x in args.ymin]
    ymax = [None if (x.lower() in ('auto', '_')) else float(x)
            for x in args.ymax]

    # Settings should be list corresponding to n_plots, or value for all
    ymin = ymin[0] if len(ymin) == 1 else ymin
    ymax = ymax[0] if len(ymax) == 1 else ymax

    optplot(modes=args.mode, filenames=args.filenames, codes=args.code,
            prefix=args.prefix, directory=args.directory,
            gaussian=args.gaussian, band_gaps=args.bandgaps,
            labels=args.labels, average=args.anisotropic, height=args.height,
            width=args.width, xmin=args.xmin, xmax=args.xmax, ymin=ymin,
            ymax=ymax, colours=None, image_format=args.image_format,
            dpi=args.dpi, style=args.style, no_base_style=args.no_base_style,
            fonts=args.font)


if __name__ == "__main__":
    main()
