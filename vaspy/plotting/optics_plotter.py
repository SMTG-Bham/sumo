# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import numpy as np

from matplotlib.ticker import MaxNLocator, FuncFormatter

from vaspy.electronic_structure.dos import sort_orbitals
from vaspy.plotting import (pretty_plot, default_colours,
                            power_tick)

line_width = 1.5
label_size = 22
optics_colours = np.array([[23, 71, 158], [217, 59, 43],
                           [13, 177, 75], [247, 148, 51],
                           [13, 177, 75]] + default_colours) / 255.


class VOpticsPlotter(object):

    def __init__(self, abs_data, band_gap=None, label=None):
        """Vaspy class for plotting optical absorption spectra.

        The class should be initialised with the absorption data
        from the optics.calculate_alpha() method.

        Args:
            abs_data (tuple or list): The optical absorption spectra,
                formatted as a tuple of ([energies], [alpha]). Alternatively,
                the ansiotropic absorption can be plotted through data
                formatted as ([energies], [alphaxx, alphayy, alphazz]).
                If a list of tuples are provided then multiple absorption
                spectra can be plotted simultaneously.
            band_gap (float or list): The fundamental band gap of the
                material to be plotted as dashed line. If plotting multiple
                spectra then a list of band gaps can be provided.
            label (str or list): A label to identify the spectra. If
                plotting multiple spectra then a list of labels can be provided.
        """
        if type(abs_data) is tuple:
            abs_data = [abs_data]

        if type(band_gap) is float:
            band_gap = [band_gap]
        elif not band_gap:
            band_gap = [None] * len(abs_data)

        if type(label) is str:
            label = [label]
        elif not label and len(abs_data) > 1:
            label = [str(i) for i in range(1, len(abs_data) + 1)]
        elif not label:
            label = [''] * len(abs_data)

        if len(set([len(label), len(abs_data), len(band_gap)])) > 1:
            raise ValueError('abs_data, band_gap, and label not same size')

        # xmax: find lowest energy where absorption > 2e4; add 1 eV
        xmax = 0
        for ener, alpha in abs_data:
            mask = alpha > 2e4
            if len(alpha.shape) == 1:
                x = min(ener[mask])
            else:
                x = max([min(ener[mask[:, i]]) for i in range(3)])
            xmax = x if x > xmax else xmax

        self._abs_data = abs_data
        self._band_gap = band_gap
        self._label = label
        self._xmax = xmax + 1.

    def get_plot(self, width=6., height=6., xmin=0., xmax=None, ymin=0,
                 ymax=1e5, colours=None, dpi=400, plt=None, fonts=None):
        """Get a matplotlib pyplot object of the density of states.

        Args:
            width (float): The width of the graph.
            height (float): The height of the graph.
            xmin (float): The minimum energy to plot.
            xmax (float): The maximum energy to plot.
            ymin (float): The minimum absorption intensity to plot.
            ymax (float): The maximum absorption intensity to plot.
            colours (dict): Specify custom colours as {'Element': colour} where
                colour is a hex number.
            dpi (int): The dots-per-inch (pixel density) for the image.
            plt (pyplot object): Matplotlib pyplot object to use for plotting.
            fonts (list): List of fonts to use in the plot.

        Returns:
            matplotlib pyplot object.
        """
        plt = pretty_plot(width=width, height=height, dpi=dpi, plt=plt,
                          fonts=fonts)
        ax = plt.gca()

        for (ener, alpha), l, bg, c in zip(self._abs_data, self._label,
                                           self._band_gap, optics_colours):
            if len(alpha.shape) == 1:
                # if averaged optics only plot one line
                ax.plot(ener, alpha, lw=line_width, label=l, c=c)
            else:
                for i, d, ls in zip(range(3), ['xx', 'yy', 'zz'],
                                    ['-', '--', '-.']):
                    n = d if l == '' else '{}$_\mathregular{{{}}}$'.format(l, d)
                    ax.plot(ener, alpha[:, i], lw=line_width, ls=ls,
                            label=n, c=c)

            # plot band gap line
            if bg:
                ax.plot([bg, bg], [ymin, ymax], lw=line_width, ls=':', c=c)

        xmax = xmax if xmax else self._xmax
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        ax.tick_params(axis='x', which='both', top='off')
        ax.tick_params(axis='x', which='both', right='off')
        ax.yaxis.set_major_formatter(FuncFormatter(power_tick))
        ax.yaxis.set_major_locator(MaxNLocator(4))

        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Absorption (cm$^\mathregular{-1}$)')

        if not np.all(np.array(self._label) == ''):
            ax.legend(loc='best', frameon=False, ncol=1,
                      prop={'size': label_size - 3})

        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))
        plt.tight_layout()
        return plt
