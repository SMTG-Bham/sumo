# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
This module provides a class for plotting phonon band structure diagrams.
"""

import itertools
import json
import logging

from matplotlib import rcParams
from matplotlib.cbook import flatten
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
from matplotlib.transforms import blended_transform_factory
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.plotter import PhononBSPlotter

from sumo.plotting import (
    draw_themed_line,
    pretty_plot,
    pretty_subplot,
    styled_plot,
    sumo_base_style,
    sumo_bs_style,
    sumo_phonon_style,
)


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

    @staticmethod
    def _plot_phonon_dos(dos, ax=None, color=None, dashline=False):
        if ax is None:
            import matplotlib.pyplot as plt

            ax = plt.gca()
        if color is None:
            color = "C0"
        y, x = dos[:, 0], dos[:, 1]
        ax.plot(x, y, "-", color=color)
        ax.fill_betweenx(y, x, 0, color=color, alpha=0.5)
        ax.set_xticks([])
        ax.set_xlim([0, max(x) * 1.1])
        ax.set_xlabel("DOS")

        if dashline:
            draw_themed_line(0, ax)

    @styled_plot(sumo_base_style, sumo_bs_style, sumo_phonon_style)
    def get_plot(
        self,
        units="THz",
        ymin=None,
        ymax=None,
        width=None,
        height=None,
        dpi=None,
        plt=None,
        fonts=None,
        dos=None,
        dos_aspect=3,
        color=None,
        style=None,
        no_base_style=False,
        from_json=None,
        legend=None,
    ):
        """Get a :obj:`matplotlib.pyplot` object of the phonon band structure.

        Args:
            units (:obj:`str`, optional): Units of phonon frequency. Accepted
                (case-insensitive) values are Thz, cm-1, eV, meV.
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
            dos (:obj:`np.ndarray`): 2D Numpy array of total DOS data
            dos_aspect (float): Width division for vertical DOS
            color (:obj:`str` or :obj:`tuple`, optional): Line/fill colour in
                any matplotlib-accepted format
            style (:obj:`list`, :obj:`str`, or :obj:`dict`): Any matplotlib
                style specifications, to be composed on top of Sumo base
                style.
            no_base_style (:obj:`bool`, optional): Prevent use of sumo base
                style. This can make alternative styles behave more
                predictably.
            from_json (:obj:`list` or :obj:`None`, optional): List of paths to
                :obj:`pymatgen.phonon.bandstructure.PhononBandStructureSymmline`
                JSON dump files. These are used to generate additional plots
                displayed under the data attached to this plotter.
                The k-point path should be the same as the main plot; the
                reciprocal lattice is adjusted to fit the scaling of the main
                data input.

        Returns:
            :obj:`matplotlib.pyplot`: The phonon band structure plot.
        """
        if from_json is None:
            from_json = []

        if legend is None:
            legend = [""] * (len(from_json) + 1)
        else:
            if len(legend) == 1 + len(from_json):
                pass
            elif len(legend) == len(from_json):
                legend = [""] + list(legend)
            else:
                raise ValueError("Inappropriate number of legend entries")

        if color is None:
            color = "C0"  # Default to first colour in matplotlib series

        if dos is not None:
            plt = pretty_subplot(
                1,
                2,
                width=width,
                height=height,
                sharex=False,
                sharey=True,
                dpi=dpi,
                plt=plt,
                gridspec_kw={"width_ratios": [dos_aspect, 1], "wspace": 0},
            )
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt)
            ax = plt.gca()

        def _plot_lines(data, ax, color=None, alpha=1, zorder=1):
            """Pull data from any PhononBSPlotter and add to axis"""
            dists = data["distances"]
            freqs = data["frequency"]

            # nd is branch index, nb is band index, nk is kpoint index
            for nd, nb in itertools.product(
                range(len(data["distances"])), range(self._bs.nb_bands)
            ):
                f = freqs[nd][nb]

                # plot band data
                ax.plot(dists[nd], f, ls="-", c=color, zorder=zorder)

        data = self.bs_plot_data()
        _plot_lines(data, ax, color=color)

        for i, bs_json in enumerate(from_json):
            with open(bs_json) as f:
                json_data = json.load(f)
                json_data["lattice_rec"] = json.loads(self._bs.lattice_rec.to_json())
                bs = PhononBandStructureSymmLine.from_dict(json_data)

                # bs.lattice_rec = self._bs.lattice_rec
                # raise Exception(bs.qpoints)
            json_plotter = PhononBSPlotter(bs)
            json_data = json_plotter.bs_plot_data()
            if json_plotter.n_bands != self._bs.nb_bands:
                raise Exception(
                    f"Number of bands in {bs_json} does not match main plot"
                )
            _plot_lines(json_data, ax, color=f"C{i + 1}", zorder=0.5)

        if any(legend):  # Don't show legend if all entries are empty string
            from matplotlib.lines import Line2D

            ax.legend(
                [Line2D([0], [0], color=f"C{i}") for i in range(len(legend))],
                legend,
            )

        self._maketicks(ax, units=units)
        self._makeplot(
            ax,
            plt.gcf(),
            data,
            width=width,
            height=height,
            ymin=ymin,
            ymax=ymax,
            dos=dos,
            color=color,
        )
        plt.tight_layout()
        plt.subplots_adjust(wspace=0)

        return plt

    def _makeplot(
        self,
        ax,
        fig,
        data,
        ymin=None,
        ymax=None,
        height=6,
        width=6,
        dos=None,
        color=None,
    ):
        """Utility method to tidy phonon band structure diagrams."""
        # Define colours
        if color is None:
            color = "C0"  # Default to first colour in matplotlib series

        # set x and y limits
        tymax = ymax if (ymax is not None) else max(flatten(data["frequency"]))
        tymin = ymin if (ymin is not None) else min(flatten(data["frequency"]))
        pad = (tymax - tymin) * 0.05

        if ymin is None:
            ymin = 0 if tymin >= self.imag_tol else tymin - pad
        ymax = ymax if ymax else tymax + pad

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(0, data["distances"][-1][-1])

        if ymin < 0:
            dashline = True
            draw_themed_line(0, ax)

        else:
            dashline = False

        if dos is not None:
            self._plot_phonon_dos(dos, ax=fig.axes[1], color=color, dashline=dashline)
        else:
            # keep correct aspect ratio; match axis to canvas
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()

            if width is None:
                width = rcParams["figure.figsize"][0]
            if height is None:
                height = rcParams["figure.figsize"][1]
            ax.set_aspect((height / width) * ((x1 - x0) / (y1 - y0)))

    def _maketicks(self, ax, units="THz"):
        """Utility method to add tick marks to a band structure."""
        # set y-ticks
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))

        # set x-ticks; only plot the unique tick labels
        ticks = self.get_ticks()
        unique_d = []
        unique_l = []
        if ticks["distance"]:
            temp_ticks = list(zip(ticks["distance"], ticks["label"]))
            unique_d.append(temp_ticks[0][0])
            unique_l.append(temp_ticks[0][1])
            for i in range(1, len(temp_ticks)):
                if unique_l[-1] != temp_ticks[i][1]:
                    unique_d.append(temp_ticks[i][0])
                    unique_l.append(temp_ticks[i][1])

        logging.info("\nLabel positions:")
        for dist, label in list(zip(unique_d, unique_l)):
            logging.info(f"\t{dist:.4f}: {label}")

        ax.set_xticks(unique_d)
        ax.set_xticklabels(unique_l)
        ax.xaxis.grid(True, ls="-")

        trans_xdata_yaxes = blended_transform_factory(ax.transData, ax.transAxes)
        ax.vlines(
            unique_d,
            0,
            1,
            transform=trans_xdata_yaxes,
            colors=rcParams["grid.color"],
            linewidth=rcParams["grid.linewidth"],
        )

        # Use a text hyphen instead of a minus sign because some nice fonts
        # like Whitney don't come with a real minus
        labels = {
            "thz": "THz",
            "cm-1": r"cm$^{\mathrm{-}\mathregular{1}}$",
            "ev": "eV",
            "mev": "meV",
        }
        ax.set_ylabel(f"Frequency ({labels[units.lower()]})")
