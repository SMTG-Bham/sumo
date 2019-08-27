# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Plot high symmetry points on the Brillouin Zone from calculated band structure

TODO:
    - Connect the high symmetry points to make a path as it appears on the band
      structure
    - Incorporate an option to open a gui to inspect bz by eye
    - Apply user styles / appearance options
    - Modify the matplotlib backend such that both figures can be saved and a
      gui can be used

"""

import os
import sys
import glob
import logging
import argparse
import warnings

from pymatgen.io.vasp.outputs import BSVasprun
from pymatgen.electronic_structure.bandstructure import get_reconstructed_band_structure
from pymatgen.electronic_structure.plotter import plot_brillouin_zone
import matplotlib as mpl
mpl.use("Agg")
from sumo.plotting import pretty_plot_3d

__author__ = "Arthur Youd"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "August 21, 2019"


def brillplot(filenames=None, prefix=None, directory=None,
              width=6, height=6, fonts=None,
              image_format="pdf", dpi=400):
    """Generate plot of first brillouin zone from a band-structure calculation.
    Args:
        filenames (:obj:`str` or :obj:`list`, optional): Path to input files.
            Vasp:
                Use vasprun.xml or vasprun.xml.gz file.
        image_format (:obj:`str`, optional): The image file format. Can be any
            format supported by matplotlib, including: png, jpg, pdf, and svg.
            Defaults to pdf.
        dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for the
            image.
    """
    if not filenames:
        filenames = find_vasprun_files()
    elif isinstance(filenames, str):
        filenames = [filenames]
    bandstructures = []
    for vr_file in filenames:
        vr = BSVasprun(vr_file)
        bs = vr.get_band_structure(line_mode=True)
        bandstructures.append(bs)
    bs = get_reconstructed_band_structure(bandstructures)

    labels = {}
    for k in bs.kpoints:
        if k.label:
            labels[k.label] = k.frac_coords

    lines = []
    for b in bs.branches:
        lines.append([bs.kpoints[b['start_index']].frac_coords,
                      bs.kpoints[b['end_index']].frac_coords])

    plt = pretty_plot_3d(width, height, dpi=dpi, fonts=fonts)
    fig = plot_brillouin_zone(bs.lattice_rec, lines=lines, labels=labels,
                              ax=plt.gca())

    basename = "brillouin.{}".format(image_format)
    filename = "{}_{}".format(prefix, basename) if prefix else basename
    if directory:
        filename = os.path.join(directory, filename)
    fig.savefig(filename, format=image_format, dpi=dpi, bbox_inches="tight")
    return plt

def find_vasprun_files():
    """Search for vasprun files from the current directory.

    The precedence order for file locations is:
    1. First search for folders named: 'split-0*'
    2. Else, look in the current directory.
    The split folder names should always be zero based, therefore easily
    sortable.

    Returns: list of str

    """
    folders = glob.glob("split-*")
    folders = sorted(folders) if folders else ["."]

    filenames = []
    for fol in folders:
        vr_file = os.path.join(fol, "vasprun.xml")
        vr_file_gz = os.path.join(fol, "vasprun.xml.gz")
        if os.path.exists(vr_file):
            filenames.append(vr_file)
        elif os.path.exists(vr_file_gz):
            filenames.append(vr_file_gz)
        else:
            logging.error("ERROR: No vasprun.xml found in {}!".format(fol))
            sys.exit()
    return filenames


def _get_parser():
    parser = argparse.ArgumentParser(description="""
    brillplot is a script to produce publication-ready
    brillouin zone diagrams""",
                                     epilog="""
    Author: {}
    Version: {}
    Last updated: {}""".format(__author__, __version__, __date__))
    parser.add_argument(
        "-f",
        "--filenames",
        default=None,
        nargs="+",
        metavar="F",
        help="one or more vasprun.xml files to plot",
    )
    parser.add_argument(
        "-d", "--directory", metavar="D", help="output directory for files"
    )
    parser.add_argument(
        "--format",
        type=str,
        default="pdf",
        dest="image_format",
        metavar="FORMAT",
        help="image file format (options: pdf, svg, jpg, png)",
    )
    parser.add_argument(
        "--dpi", type=int, default=400, help="pixel density for image file"
    )
    return parser


def main():
    args = _get_parser().parse_args()
    logging.basicConfig(
        filename="sumo-brillplot.log",
        level=logging.INFO,
        filemode="w",
        format="%(message)s",
    )
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger("").addHandler(console)
    warnings.filterwarnings("ignore",
                            category=UserWarning, module="matplotlib")
    warnings.filterwarnings("ignore",
                            category=UnicodeWarning, module="matplotlib")
    warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

    brillplot(
        filenames=args.filenames,
        directory=args.directory,
        image_format=args.image_format,
        dpi=args.dpi,
    )
