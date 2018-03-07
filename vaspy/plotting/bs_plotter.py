# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import logging
import itertools
import copy

import numpy as np

from matplotlib.ticker import MaxNLocator, AutoMinorLocator

from vaspy.plotting import pretty_plot, pretty_subplot

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

line_width = 1.5
empty_space = 1.05
label_size = 22
band_linewidth = 2


class VBSPlotter(BSPlotter):

    def __init__(self, bs):
        """Vaspy class for plotting band structures.

        This class is similar to the pymatgen BSPlotter class but overrides some
        methods to generate prettier plots.

        Further functionality, such as projected band structure plots are
        available.

        Args:
            bs (BandStructure): A pymatgen BandStructure object.
        """
        BSPlotter.__init__(self, bs)

    def get_plot(self, zero_to_efermi=True, ymin=-6., ymax=6.,
                 width=6., height=6., vbm_cbm_marker=False, dpi=400, plt=None,
                 dos_plotter=None, dos_options=None, dos_aspect=3, fonts=None):
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

        data = self.bs_plot_data(zero_to_efermi)
        dists = data['distances']
        eners = data['energy']

        if self._bs.is_spin_polarized or self._bs.is_metal():
            is_vb = True
        else:
            is_vb = self._bs.bands[Spin.up] <= self._bs.get_vbm()['energy']

        # nd is branch index, nb is band index, nk is kpoint index
        for nd, nb in itertools.product(range(len(data['distances'])),
                                        range(self._nb_bands)):
            e = eners[nd][str(Spin.up)][nb]

            # this check is very slow but works for now
            # colour valence bands blue and conduction bands orange in semiconds
            if (self._bs.is_spin_polarized or self._bs.is_metal() or
                    np.all(is_vb[nb])):
                c = '#3953A4'
            else:
                c = '#FAA316'

            # plot band data
            ax.plot(dists[nd], e, ls='-', c=c, linewidth=band_linewidth)
            if self._bs.is_spin_polarized:
                e = eners[nd][str(Spin.down)][nb]
                ax.plot(dists[nd], e, 'r--', linewidth=band_linewidth)

        self._maketicks(ax)
        self._makeplot(ax, plt.gcf(), data, zero_to_efermi=zero_to_efermi,
                       vbm_cbm_marker=vbm_cbm_marker, width=width,
                       height=height, ymin=ymin, ymax=ymax,
                       dos_plotter=dos_plotter, dos_options=dos_options)
        plt.tight_layout()
        return plt

    def get_projected_rgb_plot(self, element_orbitals, zero_to_efermi=True,
                               ymin=-6., ymax=6., width=6., height=6.,
                               vbm_cbm_marker=False, dpi=400, plt=None,
                               dos_plotter=None, dos_options=None,
                               dos_aspect=3):
        """Get a matplotlib object for a projected rgb bandstructure plot.

        The band structure line color depends on the character of the band
        (either different elemental or orbital contributions). Each element/
        orbital is associated with red, green or blue and the corresponding rgb
        color depending on the character of the band is used. The method can
        only deal with up to 3 elements/orbitals.

        Spin up and spin down are differientiated by a '-' and a '--' line

        Args:
            element_orbitals: A list of (element, orbital). For example:
                [('Bi', 's'), ('Bi', 'p'), ('S', 'p')]. If just the element is
                specified then all the orbitals of that element are summed. For
                example: [('Bi', 's'), ('Bi', 'p'), 'S']
                The order of the element/orbitals determines the colour of that
                contribution (going red, green, blue).
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
        """
        if dos_plotter:
            width = width + height/dos_aspect
            plt = pretty_subplot(1, 2, width, height, sharex=False, dpi=dpi,
                                 plt=plt,
                                 gridspec_kw={'width_ratios': [dos_aspect, 1],
                                              'wspace': 0})
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt)
            ax = plt.gca()

        if len(element_orbitals) > 3:
            raise ValueError('Too many elements/orbitals specified (max 3)')

        data = self.bs_plot_data(zero_to_efermi)
        dists = data['distances']
        eners = data['energy']
        spins = [str(Spin.up), str(Spin.down)] if \
            self._bs.is_spin_polarized else [str(Spin.up)]
        colours = self._get_ordered_projections_by_branches(element_orbitals,
                                                            data,
                                                            normalise='select')

        # ns is spin, nd is branch index, nb is band index, nk is kpoint index
        for ns, nd, nb in itertools.product(spins,
                                            range(len(data['distances'])),
                                            range(self._nb_bands)):
            colour = [c[nd][ns][nb] for c in colours]
            if len(colour) == 2:
                colour.insert(1, np.zeros((len(dists[nd]))))
            ls = '-' if ns == str(Spin.up) else '--'

            lc = self._rgbline(dists[nd], eners[nd][ns][nb], colour[0],
                               colour[1], colour[2], alpha=1,
                               linestyles=ls)
            ax.add_collection(lc)

            # alternative is to use scatter plot
            #ax.scatter(dists[nd], eners[nd][ns][nb], s=10, edgecolors='none',
            #           c=zip(colour[0], colour[1], colour[2],
            #                 np.ones(len(dists[nd]))))

        self._maketicks(ax)
        self._makeplot(ax, plt.gcf(), data, zero_to_efermi=zero_to_efermi,
                       vbm_cbm_marker=vbm_cbm_marker, ymin=ymin, ymax=ymax,
                       height=height, width=width,
                       dos_plotter=dos_plotter, dos_options=dos_options)
        # TODO: Add rgb legend
        return plt

    def get_projected_plot(self, element_orbitals, zero_to_efermi=True,
                           ymin=-6., ymax=6., width=6., height=6.,
                           vbm_cbm_marker=False, dpi=400, plt=None,
                           dos_plotter=None, dos_options=None,
                           dos_aspect=3):
        """

        Spin up and spin down are differientiated by a '-' and a '--' line

        Args:
            element_orbitals: A list of (element, orbital). For example:
                [('Bi', 's'), ('Bi', 'p'), ('S', 'p')]. If just the element is
                specified then all the orbitals of that element are summed. For
                example: [('Bi', 's'), ('Bi', 'p'), 'S']
                The order of the element/orbitals determines the colour of that
                contribution (going red, green, blue).
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
        """
        # TODO: Finish docstring
        if dos_plotter:
            width = width + height/dos_aspect
            plt = pretty_subplot(1, 2, width, height, sharex=False, dpi=dpi,
                                 plt=plt,
                                 gridspec_kw={'width_ratios': [dos_aspect, 1],
                                              'wspace': 0})
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt)
            ax = plt.gca()

        if len(element_orbitals) > 3:
            raise ValueError('Too many elements/orbitals specified (max 3)')

        data = self.bs_plot_data(zero_to_efermi)
        dists = data['distances']
        eners = data['energy']
        spins = [str(Spin.up), str(Spin.down)] if \
            self._bs.is_spin_polarized else [str(Spin.up)]
        colours = self._get_ordered_projections_by_branches(element_orbitals,
                                                            data,
                                                            normalise='select')

        # ns is spin, nd is branch index, nb is band index, nk is kpoint index
        for ns, nd, nb in itertools.product(spins,
                                            range(len(data['distances'])),
                                            range(self._nb_bands)):
            colour = [c[nd][ns][nb] for c in colours]
            if len(colour) == 2:
                colour.insert(1, np.zeros((len(dists[nd]))))
            ls = '-' if ns == str(Spin.up) else '--'

            lc = self._rgbline(dists[nd], eners[nd][ns][nb], colour[0],
                               colour[1], colour[2], alpha=1,
                               linestyles=ls)
            ax.add_collection(lc)

            # alternative is to use scatter plot
            #ax.scatter(dists[nd], eners[nd][ns][nb], s=10, edgecolors='none',
            #           c=zip(colour[0], colour[1], colour[2],
            #                 np.ones(len(dists[nd]))))

        self._maketicks(ax)
        self._makeplot(ax, plt.gcf(), data, zero_to_efermi=zero_to_efermi,
                       vbm_cbm_marker=vbm_cbm_marker, ymin=ymin, ymax=ymax,
                       height=height, width=width,
                       dos_plotter=dos_plotter, dos_options=dos_options)
        # TODO: Add rgb legend
        return plt

    def _makeplot(self, ax, fig, data, zero_to_efermi=True,
                  vbm_cbm_marker=False, ymin=-6, ymax=6, height=6, width=6,
                  dos_plotter=None, dos_options=None):
        # draw line at Fermi level if not zeroing to e-Fermi
        if not zero_to_efermi:
            ef = self._bs.efermi
            ax.axhline(ef, linewidth=2, color='k')

        # set x and y limits
        ax.set_xlim(0, data['distances'][-1][-1])
        if self._bs.is_metal() and not zero_to_efermi:
            ax.set_ylim(self._bs.efermi + ymin, self._bs.efermi + ymax)
        else:
            ax.set_ylim(ymin, ymax)

        if vbm_cbm_marker:
            for cbm in data['cbm']:
                ax.scatter(cbm[0], cbm[1], color='#D93B2B', marker='o', s=100)
            for vbm in data['vbm']:
                ax.scatter(vbm[0], vbm[1], color='#0DB14B', marker='o', s=100)

        if dos_plotter:
            ax = fig.axes[1]
            dos_options.update({'xmin': ymin, 'xmax': ymax})
            self._makedos(ax, dos_plotter, dos_options)
        else:
            # keep correct aspect ratio square
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))

    def _makedos(self, ax, dos_plotter, dos_options):
        """This is basically the same as the VDOSPlotter get_plot function."""
        plot_data = dos_plotter.dos_plot_data(**dos_options)
        mask = plot_data['mask']
        energies = plot_data['energies'][mask]
        lines = plot_data['lines']
        spins = [Spin.up] if len(lines[0][0]['dens']) == 1 else \
            [Spin.up, Spin.down]
        for i, line_set in enumerate(plot_data['lines']):
            for line, spin in itertools.product(line_set, spins):
                if spin == Spin.up:
                    label = line['label']
                    densities = line['dens'][spin][mask]
                elif spin == Spin.down:
                    label = ""
                    densities = -line['dens'][spin][mask]
                ax.fill_betweenx(energies, densities, 0, lw=0,
                                 facecolor=line['colour'],
                                 alpha=line['alpha'])
                ax.plot(densities, energies, label=label,
                        color=line['colour'], lw=line_width)

            # x and y axis reversed versus normal dos plotting
            ax.set_ylim(dos_options['xmin'], dos_options['xmax'])
            ax.set_xlim(plot_data['ymin'], plot_data['ymax'])

            ax.tick_params(axis='y', which='both', top='off')
            ax.tick_params(axis='x', which='both', labelbottom='off',
                           labeltop='off', bottom='off', top='off')

            ax.legend(loc=2, frameon=False, ncol=1,
                      prop={'size': label_size - 3},
                      bbox_to_anchor=(1., 1.))

    def _maketicks(self, ax):
        """Utility method to add tick marks to a band structure."""
        # set y-ticks
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

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
        ax.set_ylabel('Energy (eV)')

    def _get_projections_by_branches(self, dictio):
        proj = self._bs.get_projections_on_elements_and_orbitals(dictio)
        proj[Spin.up] = np.array(proj[Spin.up])
        if self._bs.is_spin_polarized:
            proj[Spin.down] = np.array(proj[Spin.down])

        proj_br = []
        for b in self._bs.branches:
            s = b['start_index']
            e = b['end_index'] + 1
            proj_br.append({str(Spin.up): proj[Spin.up][:, s:e]})
            if self._bs.is_spin_polarized:
                proj_br.append({str(Spin.down): proj[Spin.down][:, s:e]})
        return proj_br

    def _get_ordered_projections_by_branches(self, element_orbitals, data,
                                             normalise=None):
        # if we are to normalise the data later we need access to all projs
        elements = self._bs.structure.symbol_set
        orbitals = ['s', 'p', 'd', 'f']
        dictio = dict(zip(elements, orbitals*len(elements)))
        proj = self._get_projections_by_branches(dictio)

        spins = [str(Spin.up), str(Spin.down)] if \
            self._bs.is_spin_polarized else [str(Spin.up)]
        nbands = self._nb_bands

        sum_proj = []
        # now turn the proj into an list of [proj_elem1, proj_elem2, ...]
        # where proj_elem1 is array with indicies [nbranch][spin][nband][nkpt]
        proj_by_element = [[] for i in range(len(element_orbitals))]
        for nd in range(len(data['distances'])):
            nkpts = len(data['distances'][nd])
            tproj = {str(Spin.up): np.zeros((nbands, nkpts))}
            sum_p = {str(Spin.up): np.zeros((nbands, nkpts))}
            if self._bs.is_spin_polarized:
                tproj[str(Spin.down)] = np.zeros(nbands, nkpts)
                sum_p[str(Spin.down)] = np.zeros(nbands, nkpts)
            tproj = [copy.deepcopy(tproj) for i in range(len(element_orbitals))]

            for nb, nk, ns in itertools.product(range(nbands),
                                                range(nkpts), spins):
                for i, e in enumerate(element_orbitals):
                    if type(e) is tuple:
                        tproj[i][ns][nb][nk] = proj[nd][ns][nb][nk][e[0]][e[1]]
                    else:
                        tproj[i][ns][nb][nk] = np.sum([proj[nd][ns][nb][nk][e][o]
                                                       for o in ['s', 'p', 'd']],
                                                      axis=0)

                # if all then normalise against all projections
                # if select then we just normalise against specified projs
                if normalise == 'all':
                    sum_p[ns][nb][nk] = np.sum([proj[nd][ns][nb][nk][e][o]
                            for e, o in itertools.product(elements, orbitals)],
                            axis=0)
                elif normalise == 'select':
                    sum_p[ns][nb][nk] = np.sum([p[ns][nb][nk] for p in tproj])

            for i in range(len(element_orbitals)):
                if normalise:
                    # to prevent warnings/errors relating to divide by zero,
                    # catch warnings and surround divide with np.nan_to_num
                    with np.errstate(divide='ignore', invalid='ignore'):
                        for ns in spins:
                            tproj[i][ns] = np.nan_to_num(tproj[i][ns] / sum_p[ns])
                proj_by_element[i].append(tproj[i])
            sum_proj.append(sum_p)
        return proj_by_element

    @staticmethod
    def _rgbline(x, y, red, green, blue, alpha=1, linestyles="solid"):
        """An RGB colored line for plotting.

        Creation of segments based on:
        http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb

        Args:
            ax: matplotlib axis
            x: x-axis data (k-points)
            y: y-axis data (energies)
            red: red data
            green: green data
            blue: blue data
            alpha: alpha values data
            linestyles: linestyle for plot (e.g., "solid" or "dotted")
        """
        # TODO: Add interpolation
        from matplotlib.collections import LineCollection

        pts = np.array([x, y]).T.reshape(-1, 1, 2)
        seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

        nseg = len(x) - 1
        r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
        g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
        b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
        a = np.ones(nseg, np.float) * alpha
        lc = LineCollection(seg, colors=list(zip(r, g, b, a)),
                            linewidth=2, linestyles=linestyles)
        return lc
