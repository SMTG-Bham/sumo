# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
This module provides a class for plotting optical absorption spectra.
"""

import numpy as np

from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, FuncFormatter, AutoMinorLocator

from sumo.plotting import (pretty_plot, power_tick, styled_plot,
                           sumo_base_style, sumo_optics_style)


class SOpticsPlotter(object):
    """Class for plotting optical absorption spectra.

    The easiest way to initialise this class is using the output from the
    :obj:`sumo.electronic_structure.optics.calculate_alpha()` method.

    Args:
        abs_data (:obj:`tuple` or :obj:`list`): The optical absorption
            spectra. Should be formatted as a :obj:`tuple` of :obj:`list`::

                ([energies], [alpha])

            Alternatively, the anisotropic (directional dependent) absorption
            can be plotted if the data formatted as::

                ([energies], [alpha_xx, alpha_yy, alpha_zz])

            If a :obj:`list` of :obj:`tuple` is provided, then multiple
            absorption spectra can be plotted simultaneously.
        band_gap (:obj:`float` or :obj:`list`, optional): The band gap as a
            :obj:`float`, plotted as a dashed line. If plotting multiple
            spectra then a :obj:`list` of band gaps can be provided.
        label (:obj:`str` or :obj:`list`): A label to identify the spectra.
            If plotting multiple spectra then a :obj:`list` of labels can
            be provided.
    """

    def __init__(self, abs_data, band_gap=None, label=None):
        if isinstance(abs_data, tuple):
            abs_data = [abs_data]

        if isinstance(band_gap, float):
            band_gap = [band_gap]
        elif not band_gap:
            band_gap = [None] * len(abs_data)

        if isinstance(label, str):
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

    @styled_plot(sumo_base_style, sumo_optics_style)
    def get_plot(self, width=None, height=None, xmin=0., xmax=None, ymin=0,
                 ymax=1e5, colours=None, dpi=400, plt=None, fonts=None,
                 style=None, no_base_style=False):
        """Get a :obj:`matplotlib.pyplot` object of the optical spectra.

        Args:
            width (:obj:`float`, optional): The width of the plot.
            height (:obj:`float`, optional): The height of the plot.
            xmin (:obj:`float`, optional): The minimum energy on the x-axis.
            xmax (:obj:`float`, optional): The maximum energy on the x-axis.
            ymin (:obj:`float`, optional): The minimum absorption intensity on
                the y-axis.
            ymax (:obj:`float`, optional): The maximum absorption intensity on
                the y-axis.
            colours (:obj:`list`, optional): A :obj:`list` of colours to use in
                the plot. The colours can be specified as a hex code, set of
                rgb values, or any other format supported by matplotlib.
            dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
                the image.
            plt (:obj:`matplotlib.pyplot`, optional): A
                :obj:`matplotlib.pyplot` object to use for plotting.
            fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
                a single font, specified as a :obj:`str`, or several fonts,
                specified as a :obj:`list` of :obj:`str`.
            style (:obj:`list`, :obj:`str`, or :obj:`dict`): Any matplotlib
                style specifications, to be composed on top of Sumo base
                style.
            no_base_style (:obj:`bool`, optional): Prevent use of sumo base
                style. This can make alternative styles behave more
                predictably.

        Returns:
            :obj:`matplotlib.pyplot`: The plot of optical spectra.
        """
        plt = pretty_plot(width=width, height=height, dpi=dpi, plt=plt)
        ax = plt.gca()

        optics_colours = rcParams['axes.prop_cycle'].by_key()['color']
        if colours is not None:
            optics_colours = colours + optics_colours

        for (ener, alpha), abs_label, bg, c in zip(self._abs_data,
                                                   self._label,
                                                   self._band_gap,
                                                   optics_colours):
            if len(alpha.shape) == 1:
                # if averaged optics only plot one line
                ax.plot(ener, alpha, label=abs_label, c=c)

            else:
                data = zip(range(3), ['xx', 'yy', 'zz'], ['-', '--', '-.'])

                for direction_id, direction_label, ls in data:
                    if not abs_label:
                        label = direction_label
                    else:
                        label = r'{}$_\mathregular{{{}}}$'
                        label.format(direction_label, direction_id)

                    ax.plot(ener, alpha[:, direction_id], ls=ls,
                            label=label, c=c)

            if bg:
                # plot band gap line
                ax.plot([bg, bg], [ymin, ymax], ls=':', c=c)

        xmax = xmax if xmax else self._xmax
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        ax.yaxis.set_major_formatter(FuncFormatter(power_tick))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))

        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel(r'Absorption (cm$^\mathregular{-1}$)')

        if (not np.all(np.array(self._label) == '')
                or len(np.array(self._abs_data[0][1]).shape) > 1):
            ax.legend(loc='best', frameon=False, ncol=1)

        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        if width is None:
            width = rcParams['figure.figsize'][0]
        if height is None:
            height = rcParams['figure.figsize'][1]
        ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))
        plt.tight_layout()

        return plt
