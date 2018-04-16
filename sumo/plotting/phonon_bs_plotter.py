# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
This module provides a class for plotting phonon band structure diagrams.
"""

import logging
import itertools

from matplotlib.ticker import MaxNLocator
from matplotlib.cbook import flatten

from sumo.plotting import pretty_plot

from pymatgen.phonon.plotter import PhononBSPlotter

line_width = 1.5
band_linewidth = 2


class SPhononBSPlotter(PhononBSPlotter):
    """Class for plotting phonon band structures.

    This class is similar to the :obj:`pymatgen.phonon.plotter.PhononBSPlotter`
    class but overrides some methods to generate prettier plots.

    Args:
        bs (:obj:`~pymatgen.phonon.bandstructure.PhononBandStructureSymmLine`):
            The phonon band structure.
    """

    def __init__(self, bs, imag_tol=-5e-2):
        PhononBSPlotter.__init__(self, bs)
        self.imag_tol = imag_tol

    def get_plot(self, ymin=None, ymax=None, width=6., height=6., dpi=400,
                 plt=None, fonts=None):
        """Get a :obj:`matplotlib.pyplot` object of the phonon band structure.

        Args:
            ymin (:obj:`float`, optional): The minimum energy on the y-axis.
            ymax (:obj:`float`, optional): The maximum energy on the y-axis.
            width (:obj:`float`, optional): The width of the plot.
            height (:obj:`float`, optional): The height of the plot.
            dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
                the image.
            fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
                a single font, specified as a :obj:`str`, or several fonts,
                specified as a :obj:`list` of :obj:`str`.
            plt (:obj:`matplotlib.pyplot`, optional): A
                :obj:`matplotlib.pyplot` object to use for plotting.

        Returns:
            :obj:`matplotlib.pyplot`: The phonon band structure plot.
        """
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
            ax.plot(dists[nd], f, ls='-', c='#3953A4',
                    linewidth=band_linewidth)

        self._maketicks(ax)
        self._makeplot(ax, plt.gcf(), data, width=width, height=height,
                       ymin=ymin, ymax=ymax)
        plt.tight_layout()
        return plt

    def _makeplot(self, ax, fig, data, ymin=None, ymax=None, height=6,
                  width=6):
        """Utility method to tidy phonon band structure diagrams. """
        # set x and y limits
        tymax = ymax if ymax else max(flatten(data['frequency']))
        tymin = ymin if ymin else min(flatten(data['frequency']))
        pad = (tymax - tymin) * 0.05

        if not ymin:
            ymin = 0 if tymin >= self.imag_tol else tymin - pad
        ymax = ymax if ymax else tymax + pad

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(0, data['distances'][-1][-1])

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

        logging.info('\nLabel positions:')
        for dist, label in list(zip(unique_d, unique_l)):
            logging.info('\t{:.4f}: {}'.format(dist, label))

        ax.set_xticks(unique_d)
        ax.set_xticklabels(unique_l)
        ax.xaxis.grid(True, c='k', ls='-', lw=line_width)
        ax.set_ylabel('Frequency (THz)')
