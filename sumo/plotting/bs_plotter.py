# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
This module provides a class for plotting electronic band structure diagrams.
"""

import itertools as it
import logging

import numpy as np
from matplotlib import cycler, rcParams
from matplotlib.style import context
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import BSPlotter
from scipy.interpolate import interp1d

from sumo.electronic_structure.bandstructure import (
    force_branches,
    get_projections_by_branches,
)
from sumo.plotting import (
    colorline,
    draw_themed_line,
    pretty_plot,
    pretty_subplot,
    styled_plot,
    sumo_base_style,
    sumo_bs_style,
)

label_size = 22


class SBSPlotter(BSPlotter):
    """Class for plotting electronic band structures.

    This class is similar to the pymatgen
    :obj:`pymatgen.electronic_structure.plotter.BSPlotter` class but overrides
    some methods to generate prettier plots.

    Additional functionality, such as projected band structure plots are
    available.

    Args:
        bs (:obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
            The band structure.
    """

    def __init__(self, bs):
        if isinstance(bs, list):
            bs = [force_branches(b) for b in bs]
        else:
            bs = force_branches(bs)

        BSPlotter.__init__(self, bs)

        # old versions of pymatgen only support a single band structure
        if isinstance(self._bs, list):
            self.bs = self._bs[0]
            self.nbands = self._nb_bands[0]
        else:
            self.bs = self._bs
            self.nbands = self._nb_bands

    @staticmethod
    def _reset_zero_energy(bs_plot_data, zero_energy=0.0):
        """Modify reference energy of data from bs_plot_data

        This method is defined in Pymatgen and obtains plotting data with zero
        set to 0. (if zero_to_efermi=False), self.efermi (if
        zero_to_efermi=True and metallic) or the absolute VBM energy
        (zero_to_efermi=False and insulating).

        Sumo allows other reference energies to be specified (e.g. from
        sc-fermi); conveniently the data dictionary from Pymatgen contains a
        "zero_energy" field so the data can always be related to absolute (DFT)
        values.

        Args:
            bs_plot_data (:obj:`dict`): data dictionary from
                pymatgen.electronicstructure.plotter.BSPlotter.bs_plot_data()
            zero_energy (:obj:`float`, Optional): New reference energy

        Returns:
            dict in format of bs_plot_data.
        """

        shifted_data = {}
        energy_shift = bs_plot_data["zero_energy"] - zero_energy

        for key, value in bs_plot_data.items():
            if key in ("vbm", "cbm"):
                shifted_data[key] = [(pt[0], pt[1] + energy_shift) for pt in value]
            elif key == "energy":
                shifted_data["energy"] = {}
                for spin, energies in value.items():
                    shifted_data["energy"][spin] = [
                        array + energy_shift for array in energies
                    ]
            elif key == "zero_energy":
                shifted_data[key] = zero_energy
            else:
                shifted_data[key] = value
        return shifted_data

    @styled_plot(sumo_base_style, sumo_bs_style)
    def get_plot(
        self,
        zero_to_efermi=True,
        zero_line=False,
        zero_energy=None,
        ymin=-6.0,
        ymax=6.0,
        width=None,
        height=None,
        vbm_cbm_marker=False,
        ylabel="Energy (eV)",
        dpi=None,
        plt=None,
        plot_dos_legend=True,
        dos_plotter=None,
        dos_options=None,
        dos_label=None,
        dos_aspect=3,
        aspect=None,
        spin=None,
        fonts=None,
        style=None,
        no_base_style=False,
        title=None,
    ):
        """Get a :obj:`matplotlib.pyplot` object of the band structure.

        If the system is spin polarised, and no spin has been specified, orange
        lines are spin up, dashed blue lines are spin down. For metals, all
        bands are coloured blue. For semiconductors, blue lines indicate
        valence bands and orange lines indicates conduction bands.

        Args:
            zero_to_efermi (:obj:`bool`): Shift the plot such that the
                Fermi energy of the band structure data is plotted at 0 eV.
                Note that this "efermi" is for the benefit of Pymatgen and may
                not reflect the actual Fermi level; in Sumo it has usually
                already been shifted to the VBM.
            zero_line (:obj:`bool`, optional): Draw a horizontal line at zero
            zero_energy (:obj:`float`, optional): Zero energy reference. (If
                unset, defaults to VBM or Fermi Energy as appropriate.)
            ymin (:obj:`float`, optional): The minimum energy on the y-axis.
            ymax (:obj:`float`, optional): The maximum energy on the y-axis.
            width (:obj:`float`, optional): The width of the plot.
            height (:obj:`float`, optional): The height of the plot.
            vbm_cbm_marker (:obj:`bool`, optional): Plot markers to indicate
                the VBM and CBM locations.
            ylabel (:obj:`str`, optional): y-axis (i.e. energy) label/units
            dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
                the image.
            plt (:obj:`matplotlib.pyplot`, optional): A
                :obj:`matplotlib.pyplot` object to use for plotting.
            dos_plotter (:obj:`~sumo.plotting.dos_plotter.SDOSPlotter`, \
                optional): Plot the density of states alongside the band
                structure. This should be a
                :obj:`~sumo.plotting.dos_plotter.SDOSPlotter` object
                initialised with the data to plot.
            dos_options (:obj:`dict`, optional): The options for density of
                states plotting. This should be formatted as a :obj:`dict`
                containing any of the following keys:

                    "yscale" (:obj:`float`)
                        Scaling factor for the y-axis.
                    "xmin" (:obj:`float`)
                        The minimum energy to mask the energy and density of
                        states data (reduces plotting load).
                    "xmax" (:obj:`float`)
                        The maximum energy to mask the energy and density of
                        states data (reduces plotting load).
                    "colours" (:obj:`dict`)
                        Use custom colours for specific element and orbital
                        combinations. Specified as a :obj:`dict` of
                        :obj:`dict` of the colours. For example::

                            {
                                'Sn': {'s': 'r', 'p': 'b'},
                                'O': {'s': '#000000'}
                            }

                        The colour can be a hex code, series of rgb value, or
                        any other format supported by matplotlib.
                    "plot_total" (:obj:`bool`)
                        Plot the total density of states. Defaults to ``True``.
                    "legend_cutoff" (:obj:`float`)
                        The cut-off (in % of the maximum density of states
                        within the plotting range) for an elemental orbital to
                        be labelled in the legend. This prevents the legend
                        from containing labels for orbitals that have very
                        little contribution in the plotting range.
                    "subplot" (:obj:`bool`)
                        Plot the density of states for each element on separate
                        subplots. Defaults to ``False``.

            dos_label (:obj:`str`, optional): DOS axis label/units
            dos_aspect (:obj:`float`, optional): Aspect ratio for the band
                structure and density of states subplot. For example,
                ``dos_aspect = 3``, results in a ratio of 3:1, for the band
                structure:dos plots.
            plot_dos_legend (:obj:`bool`): Whether to plot the dos legend.
            aspect (:obj:`float`, optional): The aspect ratio of the band
                structure plot. By default the dimensions of the figure size
                are used to determine the aspect ratio. Set to ``1`` to force
                the plot to be square.
            spin (:obj:`Spin`, optional): Plot a spin-polarised band structure,
                "up" or "1" for spin up only, "down" or "-1" for spin down only.
                Defaults to ``None``.
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
            :obj:`matplotlib.pyplot`: The electronic band structure plot.
        """
        if dos_plotter:
            plt = pretty_subplot(
                1,
                2,
                width=width,
                height=height,
                sharex=False,
                dpi=dpi,
                plt=plt,
                gridspec_kw={"width_ratios": [dos_aspect, 1], "wspace": 0},
            )
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width=width, height=height, dpi=dpi, plt=plt)
            ax = plt.gca()

        data = self.bs_plot_data(zero_to_efermi=True)
        if zero_energy is not None:
            data = self._reset_zero_energy(data, zero_energy=zero_energy)

        dists = data["distances"]
        eners = data["energy"]

        if spin is not None and not self.bs.is_spin_polarized:
            raise ValueError(
                "Spin-selection only possible with spin-polarised "
                "calculation results"
            )
        elif self.bs.is_metal() or (self.bs.is_spin_polarized and not spin):
            # if metal or spin polarized and spin not specified
            is_vb = [True]
        elif spin:
            # not metal, spin-polarized and spin is set
            is_vb = self.bs.bands[spin] <= self.bs.get_vbm()["energy"]
        else:
            # not metal, not spin polarized and therefore spin not set
            is_vb = self.bs.bands[Spin.up] <= self.bs.get_vbm()["energy"]

        # nd is branch index, nb is band index, nk is kpoint index
        for nd, nb in it.product(range(len(data["distances"])), range(self.nbands)):
            e = (
                eners[str(spin)][nd][nb]
                if spin is not None
                else eners[str(Spin.up)][nd][nb]
            )

            # For closed-shell calculations with a bandgap, colour valence
            # bands blue (C0) and conduction bands orange (C1)
            #
            # For closed-shell calculations with no bandgap, colour with C0
            #
            # For spin-polarized calculations, colour spin up channel with C1
            # and overlay with C0 (dashed) spin down channel

            if self.bs.is_spin_polarized and spin is None:
                c = "C1"
            elif self.bs.is_metal() or np.all(is_vb[nb]):
                c = "C0"
            else:
                c = "C1"

            ax.plot(dists[nd], e, ls="-", c=c, zorder=1)

        # Plot second spin channel if it exists and no spin selected
        if self.bs.is_spin_polarized and spin is None:
            for nd, nb in it.product(range(len(data["distances"])), range(self.nbands)):
                e = eners[str(Spin.down)][nd][nb]
                ax.plot(dists[nd], e, c="C0", linestyle="--", zorder=2)

        self._maketicks(ax, ylabel=ylabel)

        self._makeplot(
            ax,
            plt.gcf(),
            data,
            zero_to_efermi=zero_to_efermi,
            zero_line=zero_line,
            vbm_cbm_marker=vbm_cbm_marker,
            width=width,
            height=height,
            ymin=ymin,
            ymax=ymax,
            dos_plotter=dos_plotter,
            dos_options=dos_options,
            plot_dos_legend=plot_dos_legend,
            dos_label=dos_label,
            aspect=aspect,
            spin=spin,
            title=title,
        )
        return plt

    @styled_plot(sumo_base_style, sumo_bs_style)
    def get_projected_plot(
        self,
        selection,
        mode="rgb",
        normalise="all",
        interpolate_factor=4,
        color1="#FF0000",
        color2="#0000FF",
        color3="#00FF00",
        colorspace="lab",
        circle_size=150,
        projection_cutoff=0.001,
        zero_energy=None,
        zero_to_efermi=True,
        zero_line=False,
        ymin=-6.0,
        ymax=6.0,
        width=None,
        height=None,
        vbm_cbm_marker=False,
        ylabel="Energy (eV)",
        dpi=400,
        plt=None,
        dos_plotter=None,
        dos_options=None,
        dos_label=None,
        plot_dos_legend=True,
        dos_aspect=3,
        aspect=None,
        fonts=None,
        style=None,
        no_base_style=False,
        spin=None,
        title=None,
    ):
        """Get a :obj:`matplotlib.pyplot` of the projected band structure.

        If the system is spin polarised, no spin has been specified and
        ``mode = 'rgb'`` spin up and spin down bands are differentiated by
        solid and dashed lines, respectively.
        For the other modes, spin up and spin down are plotted separately.

        Args:
            selection (list): A list of :obj:`tuple` or :obj:`string`
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

            normalise (:obj:`str`, optional): Normalisation the projections.
                Options are:

                  * ``'all'``: Projections normalised against the sum of all
                       other projections.
                  * ``'select'``: Projections normalised against the sum of the
                       selected projections.
                  * ``None``: No normalisation performed.

            interpolate_factor (:obj:`int`, optional): The factor by which to
                interpolate the band structure (necessary to make smooth
                lines). A larger number indicates greater interpolation.
            circle_size (:obj:`float`, optional): The area of the circles used
                when ``mode = 'stacked'``.
            color1 (str): A color specified in any way supported by matplotlib. Used
                when ``mode = 'rgb'``.
            color2 (str): A color specified in any way supported by matplotlib. Used
                when ``mode = 'rgb'``.
            color3 (str): A color specified in any way supported by matplotlib. Used
                when ``mode = 'rgb'``.
            colorspace (str): The colorspace in which to perform the interpolation. The
                allowed values are rgb, hsv, lab, luvlc, lablch, and xyz. Used
                when ``mode = 'rgb'``.
            projection_cutoff (:obj:`float`): Don't plot projections with
                intensities below this number. This option is useful for
                stacked plots, where small projections clutter the plot.
            zero_to_efermi (:obj:`bool`): Normalise the plot such that the
                valence band maximum is set as 0 eV.
            zero_line (:obj:`bool`, optional): Draw a horizontal line at zero
            zero_energy (:obj:`float`, optional): Zero energy reference. (If
                unset, defaults to VBM or Fermi Energy as appropriate.)
            ymin (:obj:`float`, optional): The minimum energy on the y-axis.
            ymax (:obj:`float`, optional): The maximum energy on the y-axis.
            width (:obj:`float`, optional): The width of the plot.
            height (:obj:`float`, optional): The height of the plot.
            vbm_cbm_marker (:obj:`bool`, optional): Plot markers to indicate
                the VBM and CBM locations.
            ylabel (:obj:`str`, optional): y-axis (i.e. energy) label/units
            dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
                the image.
            plt (:obj:`matplotlib.pyplot`, optional): A
                :obj:`matplotlib.pyplot` object to use for plotting.
            dos_plotter (:obj:`~sumo.plotting.dos_plotter.SDOSPlotter`, \
                optional): Plot the density of states alongside the band
                structure. This should be a
                :obj:`~sumo.plotting.dos_plotter.SDOSPlotter` object
                initialised with the data to plot.
            dos_options (:obj:`dict`, optional): The options for density of
                states plotting. This should be formatted as a :obj:`dict`
                containing any of the following keys:

                    "yscale" (:obj:`float`)
                        Scaling factor for the y-axis.
                    "xmin" (:obj:`float`)
                        The minimum energy to mask the energy and density of
                        states data (reduces plotting load).
                    "xmax" (:obj:`float`)
                        The maximum energy to mask the energy and density of
                        states data (reduces plotting load).
                    "colours" (:obj:`dict`)
                        Use custom colours for specific element and orbital
                        combinations. Specified as a :obj:`dict` of
                        :obj:`dict` of the colours. For example::

                           {
                                'Sn': {'s': 'r', 'p': 'b'},
                                'O': {'s': '#000000'}
                            }

                        The colour can be a hex code, series of rgb value, or
                        any other format supported by matplotlib.
                    "plot_total" (:obj:`bool`)
                        Plot the total density of states. Defaults to ``True``.
                    "legend_cutoff" (:obj:`float`)
                        The cut-off (in % of the maximum density of states
                        within the plotting range) for an elemental orbital to
                        be labelled in the legend. This prevents the legend
                        from containing labels for orbitals that have very
                        little contribution in the plotting range.
                    "subplot" (:obj:`bool`)
                        Plot the density of states for each element on separate
                        subplots. Defaults to ``False``.

            dos_label (:obj:`str`, optional): DOS axis label/units
            plot_dos_legend (:obj:`bool`): Whether to plot the dos legend.
            dos_aspect (:obj:`float`, optional): Aspect ratio for the band
                structure and density of states subplot. For example,
                ``dos_aspect = 3``, results in a ratio of 3:1, for the band
                structure:dos plots.
            aspect (:obj:`float`, optional): The aspect ratio of the band
                structure plot. By default the dimensions of the figure size
                are used to determine the aspect ratio. Set to ``1`` to force
                the plot to be square.
            fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
                a single font, specified as a :obj:`str`, or several fonts,
                specified as a :obj:`list` of :obj:`str`.
            style (:obj:`list`, :obj:`str`, or :obj:`dict`): Any matplotlib
                style specifications, to be composed on top of Sumo base
                style.
            no_base_style (:obj:`bool`, optional): Prevent use of sumo base
                style. This can make alternative styles behave more
                predictably.
            spin (:obj:`Spin`, optional): Plot a spin-polarised band structure,
                "up" or "1" for spin up only, "down" or "-1" for spin down
                only. Defaults to ``None``.

        Returns:
            :obj:`matplotlib.pyplot`: The projected electronic band structure
            plot.
        """
        if mode == "rgb" and len(selection) > 3:
            raise ValueError("Too many elements/orbitals specified (max 3)")
        elif mode == "solo" and dos_plotter:
            raise ValueError("Solo mode plotting with DOS not supported")

        if dos_plotter:
            plt = pretty_subplot(
                1,
                2,
                width,
                height,
                sharex=False,
                dpi=dpi,
                plt=plt,
                gridspec_kw={"width_ratios": [dos_aspect, 1], "wspace": 0},
            )
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt)
            ax = plt.gca()

        data = self.bs_plot_data(zero_to_efermi=zero_to_efermi)
        if zero_energy is not None:
            data = self._reset_zero_energy(data, zero_energy=zero_energy)

        nbranches = len(data["distances"])

        # Ensure we do spin up first, then spin down
        spins = sorted(self.bs.bands.keys(), key=lambda s: -s.value)
        if spin is not None and len(spins) == 1:
            raise ValueError(
                "Spin-selection only possible with spin-polarised "
                "calculation results"
            )

        if spin is Spin.up:
            spins = [spins[0]]
        elif spin is Spin.down:
            spins = [spins[1]]

        proj = get_projections_by_branches(self.bs, selection, normalise=normalise)

        # nd is branch index
        for spin, nd in it.product(spins, range(nbranches)):
            # mask data to reduce plotting load
            bands = np.array(data["energy"][str(spin)][nd])
            mask = np.where(
                np.any(bands > ymin - 0.05, axis=1)
                & np.any(bands < ymax + 0.05, axis=1)
            )
            distances = data["distances"][nd]
            bands = bands[mask]
            weights = [proj[nd][i][spin][mask] for i in range(len(selection))]

            if len(distances) > 2:  # Only interpolate if it makes sense to do so
                # interpolate band structure to improve smoothness
                temp_dists = np.linspace(
                    distances[0], distances[-1], len(distances) * interpolate_factor
                )
                bands = interp1d(
                    distances,
                    bands,
                    axis=1,
                    bounds_error=False,
                    fill_value="extrapolate",
                )(temp_dists)
                weights = interp1d(
                    distances,
                    weights,
                    axis=2,
                    bounds_error=False,
                    fill_value="extrapolate",
                )(temp_dists)
                distances = temp_dists

            else:  # change from list to array if we skipped the scipy interpolation
                weights = np.array(weights)
                bands = np.array(bands)
                distances = np.array(distances)

            # sometimes VASP produces very small negative weights
            weights[weights < 0] = 0

            if mode == "rgb":
                # colours aren't used now but needed later for legend
                colours = [color1, color2, color3]

                # if only two orbitals then just use red and blue
                if len(weights) == 2:
                    weights = np.insert(weights, 2, np.zeros(weights[0].shape), axis=0)
                    colours = [color1, color2]

                ls = "-" if spin == Spin.up else "--"
                lc = colorline(
                    distances,
                    bands,
                    weights.transpose((1, 2, 0)),
                    color1=color1,
                    color2=color2,
                    color3=color3,
                    colorspace=colorspace,
                    linestyles=ls,
                    linewidth=(rcParams["lines.linewidth"] * 1.25),
                )
                ax.add_collection(lc)

            elif mode == "stacked":
                # TODO: Handle spin

                # use some nice custom colours first, then default colours
                colours = ["#3952A3", "#FAA41A", "#67BC47", "#6ECCDD", "#ED2025"]
                colour_series = rcParams["axes.prop_cycle"].by_key()["color"]
                colours.extend(colour_series)

                # very small circles look crap
                weights[weights < projection_cutoff] = 0

                distances = list(distances) * len(bands)
                bands = bands.flatten()
                zorders = range(-len(weights), 0)
                for w, c, z in zip(weights, colours, zorders):
                    ax.scatter(
                        distances,
                        bands,
                        c=c,
                        s=circle_size * w**2,
                        zorder=z,
                        rasterized=True,
                    )

        # plot the legend
        for c, spec in zip(colours, selection):
            if isinstance(spec, str):
                label = spec
            else:
                label = f"{spec[0]} ({' + '.join(spec[1])})"
            ax.scatter([-10000], [-10000], c=c, s=50, label=label, edgecolors="none")

        if dos_plotter:
            loc = 1
            anchor_point = (-0.2, 1)
        else:
            loc = 2
            anchor_point = (0.95, 1)

        ax.legend(
            bbox_to_anchor=anchor_point,
            loc=loc,
            frameon=False,
            handletextpad=0.1,
            scatterpoints=1,
        )

        # finish and tidy plot
        self._maketicks(ax, ylabel=ylabel)
        self._makeplot(
            ax,
            plt.gcf(),
            data,
            zero_to_efermi=zero_to_efermi,
            zero_line=zero_line,
            vbm_cbm_marker=vbm_cbm_marker,
            width=width,
            height=height,
            ymin=ymin,
            ymax=ymax,
            dos_plotter=dos_plotter,
            dos_options=dos_options,
            dos_label=dos_label,
            plot_dos_legend=plot_dos_legend,
            aspect=aspect,
            spin=spin,
            title=title,
        )
        return plt

    def _makeplot(
        self,
        ax,
        fig,
        data,
        zero_to_efermi=True,
        zero_line=False,
        vbm_cbm_marker=False,
        ymin=-6.0,
        ymax=6.0,
        height=None,
        width=None,
        dos_plotter=None,
        dos_options=None,
        dos_label=None,
        plot_dos_legend=True,
        aspect=None,
        spin=None,
        title=None,
    ):
        """Tidy the band structure & add the density of states if required."""
        if zero_line:
            draw_themed_line(0, ax)

        # Draw line at Fermi level if not zeroing to e-Fermi
        elif not zero_to_efermi:
            ef = self.bs.efermi
            draw_themed_line(ef, ax)

        # set x and y limits
        ax.set_xlim(0, data["distances"][-1][-1])
        if self.bs.is_metal() and not zero_to_efermi:
            ax.set_ylim(self.bs.efermi + ymin, self.bs.efermi + ymax)
        else:
            ax.set_ylim(ymin, ymax)

        if vbm_cbm_marker:
            for cbm in data["cbm"]:
                ax.scatter(cbm[0], cbm[1], color="C2", marker="o")
            for vbm in data["vbm"]:
                ax.scatter(vbm[0], vbm[1], color="C3", marker="o")

        if title:
            fig.suptitle(title)

        if dos_plotter:
            ax = fig.axes[1]

            if not dos_options:
                dos_options = {}

            dos_options.update(
                {
                    "xmin": ymin,
                    "xmax": ymax,
                    "zero_energy": data["zero_energy"],
                    "zero_to_efermi": False,
                }
            )
            self._makedos(
                ax,
                dos_plotter,
                dos_options,
                dos_label=dos_label,
                plot_legend=plot_dos_legend,
                zero_line=zero_line,
                spin=spin,
            )
        else:
            # keep correct aspect ratio for axes based on canvas size
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            if width is None:
                width = rcParams["figure.figsize"][0]
            if height is None:
                height = rcParams["figure.figsize"][1]

            if aspect is not False:
                if not aspect:
                    aspect = height / width

                ax.set_aspect(aspect * ((x1 - x0) / (y1 - y0)))

    @staticmethod
    def _makedos(
        ax,
        dos_plotter,
        dos_options,
        dos_label=None,
        plot_legend=True,
        zero_line=False,
        spin=None,
    ):
        """This is basically the same as the SDOSPlotter get_plot function."""

        # don't use first 4 colours; these are the band structure line colours
        cycle = cycler("color", rcParams["axes.prop_cycle"].by_key()["color"][4:])
        with context({"axes.prop_cycle": cycle}):
            plot_data = dos_plotter.dos_plot_data(**dos_options)

        mask = plot_data["mask"]
        energies = plot_data["energies"][mask]
        lines = plot_data["lines"]
        if spin is None:
            spins = [Spin.up] if len(lines[0][0]["dens"]) == 1 else [Spin.up, Spin.down]
        elif isinstance(spin, Spin):
            spins = [spin]
        else:
            spins = spin

        # disable y ticks for DOS panel
        ax.tick_params(axis="y", which="both", right=False)

        for line_set in plot_data["lines"]:
            for line, spin in it.product(line_set, spins):
                if spin == Spin.up or len(spins) == 1:
                    label = line["label"]
                    densities = line["dens"][spin][mask]
                else:
                    label = ""
                    densities = -line["dens"][spin][mask]
                ax.fill_betweenx(
                    energies,
                    densities,
                    0,
                    lw=0,
                    facecolor=line["colour"],
                    alpha=line["alpha"],
                )
                ax.plot(densities, energies, label=label, color=line["colour"])

        # x and y axis reversed versus normal dos plotting
        ax.set_ylim(dos_options["xmin"], dos_options["xmax"])
        if len(spins) == 1:
            ax.set_xlim(0, plot_data["ymax"])
        else:
            ax.set_xlim(plot_data["ymin"], plot_data["ymax"])

        if zero_line:
            draw_themed_line(0, ax)

        if dos_label is not None:
            ax.set_xlabel(dos_label)

        ax.set_xticklabels([])
        if plot_legend:
            ax.legend(loc=2, frameon=False, ncol=1, bbox_to_anchor=(1.0, 1.0))

    @staticmethod
    def _sanitise_label(label):
        """Implement label hacks: Hide trailing @, remove label with leading @"""

        import re

        if re.match("^@.*$", label):
            return None
        else:
            return re.sub("@+$", "", label)

    @classmethod
    def _sanitise_label_group(cls, labelgroup):
        r"""Implement label hacks: Hide trailing @, remove label with leading @

        Labels split with $\mid$ symbol will be treated for each part.
        """

        if r"$\mid$" in labelgroup:
            label_components = labelgroup.split(r"$\mid$")
            good_labels = [
                i for i in map(cls._sanitise_label, label_components) if i is not None
            ]
            if len(good_labels) == 0:
                return None
            else:
                return r"$\mid$".join(good_labels)
        else:
            return cls._sanitise_label(labelgroup)

    def _maketicks(self, ax, ylabel="Energy (eV)"):
        """Utility method to add tick marks to a band structure."""
        # set y-ticks
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

        # set x-ticks; only plot the unique tick labels
        ticks = self.get_ticks()
        unique_d = []
        unique_l = []
        if ticks["distance"]:
            temp_ticks = list(zip(ticks["distance"], ticks["label"]))
            unique_d.append(temp_ticks[0][0])
            unique_l.append(temp_ticks[0][1])
            for i in range(1, len(temp_ticks)):
                # Hide labels marked with @
                if "@" in temp_ticks[i][1]:
                    # If a branch connection, check all parts of label
                    if r"$\mid$" in temp_ticks[i][1]:
                        label_components = temp_ticks[i][1].split(r"$\mid$")
                        good_labels = [i for i in label_components if i[0] != "@"]
                        if len(good_labels) == 0:
                            continue
                        else:
                            temp_ticks[i] = (
                                temp_ticks[i][0],
                                r"$\mid$".join(good_labels),
                            )
                    # If a single label, check first character
                    elif temp_ticks[i][1][0] == "@":
                        continue

                # Append label to sequence if it is not same as predecessor
                if unique_l[-1] != temp_ticks[i][1]:
                    unique_d.append(temp_ticks[i][0])
                    unique_l.append(temp_ticks[i][1])

        logging.info("Label positions:")
        for dist, label in list(zip(unique_d, unique_l)):
            logging.info(f"\t{dist:.4f}: {label}")

        ax.set_xticks(unique_d)
        ax.set_xticklabels(unique_l)
        ax.xaxis.grid(True)
        ax.set_ylabel(ylabel)
