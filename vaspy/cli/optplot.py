# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Script to calculate and plot optical spectra from a VASP calculation.
"""

__author__ = "Alex Ganose"
__version__ = "0.2"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Jan 10, 2018"

import argparse

import matplotlib as mpl
mpl.use('Agg')

from pymatgen.io.vasp import Vasprun

from vaspy.electronic_structure.plotter import VOpticsPlotter
from vaspy.electronic_structure.optics import broaden_eps, calculate_alpha, write_files


def optplot(filenames='vasprun.xml', prefix=None, directory=None,
            gaussian=None, band_gaps=None, labels=None, average=True, height=6,
            width=6, xmin=0, xmax=None, ymin=0, ymax=1e5, colours=None,
            image_format='pdf', dpi=400, plt=None, fonts=None):
    """A script to plot optical absorption spectra from VASP calculations.

    Args:
        filenames (str or list): A vasprun.xml file to plot (can be gziped).
            Alternatively, a list of vasprun.xml files can be provided, in
            which case the absorption spectra for each will be plotted.
        prefix (str): A prefix for the files generated.
        directory (str): Specify a directory in which the files are saved.
        gaussian (float): The sigma of the Gaussian broadening to apply.
        band_gaps (float or list): The fundamental band gap of the
            material to be plotted as dashed line. If plotting multiple
            spectra then a list of band gaps can be provided.
        label (str or list): A label to identify the spectra. If
            plotting multiple spectra then a list of labels can be provided.
        average (bool, optional): Average the dielectric response across
            all lattice directions.
        height (float): The height of the graph.
        width (float): The width of the graph.
        xmin (float): The minimum energy to plot.
        xmax (float): The maximum energy to plot.
        ymin (float): The minimum absorption intensity to plot.
        ymax (float): The maximum absorption intensity to plot.
        colours (dict): Specify custom colours - currently not implemented.
        image_format (str): The image file format (matplotlib only). Can be
            any format supported by matplot, including: png, jpg, pdf, and svg.
        dpi (int): The dots-per-inch (pixel density) for the image.
        plt (pyplot object): Matplotlib pyplot object to use for plotting.
        fonts (list): List of fonts to use in the plot.

    Returns:
        A matplotlib pyplot object.
    """

    if type(filenames) is str:
        filenames = [filenames]

    vrs = [Vasprun(f) for f in filenames]
    dielectrics = [vr.dielectric for vr in vrs]

    if gaussian:
        dielectrics = [broaden_eps(d, gaussian)
                       for d in dielectrics]

    abs_data = [calculate_alpha(d, average=average)
                for d in dielectrics]

    if type(band_gaps) is list and not band_gaps:
        # empty list therefore get bandgap from vasprun files
        band_gaps = [vr.get_band_structure().get_band_gap()['energy']
                     for vr in vrs]
    elif type(band_gaps) is list and 'vasprun' in band_gaps[0]:
        # band_gaps contains list of vasprun files
        bg_vrs = [Vasprun(f) for f in band_gaps]
        band_gaps = [vr.get_band_structure().get_band_gap()['energy']
                     for vr in bg_vrs]
    elif type(band_gaps) is list:
        # band_gaps is non empty list w. no vaspruns; presume floats
        band_gaps = [float(i) for i in band_gaps]

    save_files = False if plt else True

    plotter = VOpticsPlotter(abs_data, band_gap=band_gaps, label=labels)
    plt = plotter.get_plot(width=width, height=height, xmin=xmin,
                           xmax=xmax, ymin=ymin, ymax=ymax,
                           colours=colours, dpi=dpi, plt=plt, fonts=fonts)

    if save_files:
        basename = 'absorption.{}'.format(image_format)
        filename = '{}_{}'.format(prefix, basename) if prefix else basename
        if directory:
            filename = os.path.join(directory, filename)
        plt.savefig(filename, format=image_format, dpi=dpi)
        write_files(abs_data, prefix=prefix, directory=directory)
    else:
        return plt


def main():
    parser = argparse.ArgumentParser(description="""
    optics is a script to help calculate and plot optical absorption spectra
    from VASP calculations. Absorption given in units of cm^-1.""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))

    parser.add_argument('-f', '--filenames', help="""vasprun.xml file to plot.
                        Can supply more than one vasprun file to process and
                        plot multiple spectra simultaneously.""",
                        default='vasprun.xml', nargs='+')
    parser.add_argument('-p', '--prefix', help='prefix for the files generated')
    parser.add_argument('-d', '--directory', help='output directory for files')
    parser.add_argument('-g', '--gaussian', type=float,
                        help='Amount of gaussian broadening to apply')
    parser.add_argument('-b', '--bandgaps', nargs='*',
                        help="""Specify fundamental band gaps, if this option is called with
                        no arguments, the band gap for each system will be read from the
                        vasprun.xml file. Alternatively, the path to a vasprun.xml file can
                        be specified, in which case the band gap will be read from this output.
                        Furthermore, a number can be provided, in which case this value will be
                        used. If plotting multiple optical spectra then an equivalent number of
                        band gaps should be specified for this option also.""")
    parser.add_argument('-l', '--labels', nargs='+',
                        help='Labels for the absorption specta.')
    parser.add_argument('-a', '--anisotropic', action='store_false',
                        help='Give the absorption separated into to the x, y, and z directions')
    parser.add_argument('--height', type=float, default=6.,
                        help='The height of the graph')
    parser.add_argument('--width', type=float, default=6.,
                        help='The width of the graph')
    parser.add_argument('--xmin', type=float, default=0.,
                        help='The minimum energy on the x axis')
    parser.add_argument('--xmax', type=float, default=None,
                        help='The maximum energy on the x axis')
    parser.add_argument('--ymin', type=float, default=0.,
                        help='The minimum energy on the y axis')
    parser.add_argument('--ymax', type=float, default=1e5,
                        help='The maximum energy on the y axis')
    parser.add_argument('--format', type=str, default='pdf',
                        dest='image_format',
                        help='select image format from pdf, svg, jpg, & png')
    parser.add_argument('--dpi', type=int, default=400,
                        help='pixel density for generated images')
    parser.add_argument('--font', default=None, help='Font to use.')
    args = parser.parse_args()

    optplot(filenames=args.filenames, prefix=args.prefix, directory=args.directory,
            gaussian=args.gaussian, band_gaps=args.bandgaps, labels=args.labels,
            average=args.anisotropic, height=args.height, width=args.width,
            xmin=args.xmin, xmax=args.xmax, ymin=args.ymin, ymax=args.ymax,
            colours=None, image_format=args.image_format, dpi=args.dpi, fonts=[args.font])


if __name__ == "__main__":
    main()
