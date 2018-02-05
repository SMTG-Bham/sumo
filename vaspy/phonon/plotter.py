# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import logging
import itertools
import copy

import numpy as np

from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.cbook import flatten

from vaspy.electronic_structure.dos import sort_orbitals
from vaspy.misc.plotting import (pretty_plot, pretty_subplot,
                                 colour_cycle, default_colours,
                                 power_tick)

from pymatgen.phonon.plotter import PhononBSPlotter

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

line_width = 1.5
empty_space = 1.05
label_size = 22
band_linewidth = 2
col_cycle = colour_cycle()
optics_colours = np.array([[23, 71, 158], [217, 59, 43],
                           [13, 177, 75], [247, 148, 51],
                           [13, 177, 75]] + default_colours) / 255.

class VPhononBSPlotter(PhononBSPlotter):

    def __init__(self, bs, imag_tol=-5e-2):
        """Vaspy class for plotting phonon band structures.

        This class is similar to the pymatgen PhononBSPlotter class but
        overrides some methods to generate prettier plots.

        Args:
            bs (BandStructure): A pymatgen PhononBandStructureSymmLine object.
        """
        PhononBSPlotter.__init__(self, bs)
        self.imag_tol = imag_tol

    def get_plot(self, ymin=None, ymax=None, width=6., height=6., dpi=400,
                 plt=None, dos_plotter=None, dos_options=None, dos_aspect=3, 
                 fonts=None):
        """Get a matplotlib object for the bandstructure plot.

        If spin polarised or metallic, blue lines are spin up, red lines are
        spin down. Otherwise, for semiconductors, blue lines indicate valence
        bands and orange indicates conduction bands.

        Args:
            zero_to_efermi (bool): Automatically subtract off the Fermi energy
                from the eigenvalues and plot (E-Ef).
            ymin (float): The y-axis (energy) minimum limit.
            ymax (float): The y-axis (energy) maximum limit.
            width (float): The width of the figure.
            height (float): The height of the figure.
            vbm_cbm_marker (bool): Plot markers to indicate the VBM and CBM
                locations.
            dpi (int): The dots-per-inch (pixel density) for the image.
            plt (pyplot object): Matplotlib pyplot object to use for plotting.
            fonts (list): List of fonts to use in the plot.
        """
        if dos_plotter:
            width = width + height/dos_aspect
            plt = pretty_subplot(1, 2, width, height, sharex=False, dpi=dpi,
                                 plt=plt, fonts=fonts,
                                 gridspec_kw={'width_ratios': [dos_aspect, 1],
                                              'wspace': 0})
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt, fonts=fonts)
            ax = plt.gca()

        data = self.bs_plot_data()
        dists = data['distances']
        freqs = data['frequency']

        # nd is branch index, nb is band index, nk is kpoint index
        for nd, nb in itertools.product(range(len(data['distances'])),
                                        range(self._nb_bands)):
            f = freqs[nd][nb]

            # plot band data
            ax.plot(dists[nd], f, ls='-', c='#3953A4', linewidth=band_linewidth)

        self._maketicks(ax)
        self._makeplot(ax, plt.gcf(), data, width=width, height=height,
                       ymin=ymin, ymax=ymax, dos_plotter=dos_plotter,
                       dos_options=dos_options)
        plt.tight_layout()
        return plt

    def _makeplot(self, ax, fig, data, ymin=None, ymax=None, height=6, width=6,
                  dos_plotter=None, dos_options=None):
        # set x and y limits
        tymax = ymax if ymax else max(flatten(data['frequency']))
        tymin = ymin if ymin else min(flatten(data['frequency']))
        pad = (tymax - tymin) * 0.05

        if not ymin:
            ymin = 0 if tymin >= self.imag_tol else tymin - pad
        ymax = ymax if ymax else tymax + pad

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(0, data['distances'][-1][-1])

        if dos_plotter:
            ax = fig.axes[1]
            dos_options.update({'xmin': ymin, 'xmax': ymax})
            self._makedos(ax, dos_plotter, dos_options)
        else:
            # keep correct aspect ratio square
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))

    def _maketicks(self, ax):
        """Utility method to add tick marks to a band structure."""
        # set y-ticks
        ax.yaxis.set_major_locator(MaxNLocator(6))

        # set x-ticks; only plot the unique tick labels
        ticks = self.get_ticks()
        unique_d = []
        unique_l = []
        if ticks['distance']:
            temp_ticks = list(zip(ticks['distance'], ticks['label']))
            unique_d.append(temp_ticks[0][0])
            unique_l.append(temp_ticks[0][1])
            for i in range(1, len(temp_ticks)):
                if unique_l[-1] != temp_ticks[i][1]:
                    unique_d.append(temp_ticks[i][0])
                    unique_l.append(temp_ticks[i][1])

        logging.info('Label positions:')
        for dist, label in list(zip(unique_d, unique_l)):
            logging.info('\t{:.4f}: {}'.format(dist, label))

        ax.set_xticks(unique_d)
        ax.set_xticklabels(unique_l)
        ax.xaxis.grid(True, c='k', ls='-', lw=line_width)
        ax.set_ylabel('Frequency (THz)')
