import logging
import itertools
import vaspy.misc.plotting

from matplotlib.ticker import MultipleLocator

from vaspy.electronic_structure.dos import sort_orbitals
from vaspy.misc.plotting import pretty_plot, pretty_subplot, colour_cycle

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

line_width = 1.5
empty_space = 1.05
label_size = 22
col_cycle = colour_cycle()


class VDOSPlotter(object):
    """Vaspy class for plotting DOSs.

    The class should be initialised with the total DOS and partial density of
    states. The PDOS is usually generated as:

        pdos = vaspy.electronic_structure.dos.get_pdos()

    Args:
        dos (Dos): A Dos object containing the total density of states.
        pdos (dict): A dict mapping the elements and their orbitals to plot
            to Dos objects. For example:
            {'Bi': {'s': Dos, 'p': Dos}, 'S': {'s' Dos, ...}.

            Usually generated ysing the dos.get_pdos() function.
    """
    def __init__(self, dos, pdos=None):
        self._dos = dos
        self._pdos = pdos

    def dos_plot_data(self, yscale=1, xmin=-6, xmax=6, colours=None,
                      plot_total=True, legend_cutoff=3, subplot=False):
        """Plot the density of states either using matplotlib or xmgrace.

        Args:
            yscale (dict): Scaling factor for the y-axis.
            xmin (float): The minimum energy for determining energy mask.
            xmax (float): The maximum energy for determining energy mask.
            colours (dict): Specify custom colours as {'Element': colour} where
                colour is a hex number.
            plot_total (bool): Whether or not to plot total DOS.
            legend_cutoff (int): The cut-off (in % of maximum DOS plotted) for a
                elemental/orbital DOS label to appear in the legend.
            subplot (bool): Split the plot up into separate plots for each
                element.

        Returns:
            A dict with the following keys:
                energies (numpy.Array): The energies.
                mask: A numpy mask which is used to trim the density arrays and
                    prevent unwanted data being included in the output file.
                lines (dict): A list of dictionaries containing information
                    about the densities to plot. Each line contains the keys:
                        label (str): The label for the legend.
                        dens (numpy.Array): The density of states.
                        colour (str): The colour of the line.
                        alpha (float): The fill alpha.
                ymin: The minimum y-axis limit.
                ymin: The maximum y-axis limit.
        """
        # mask needed to prevent unwanted data in pdf and for finding y limit
        dos = self._dos
        pdos = self._pdos
        mask = (dos.energies >= xmin - 0.05) & (dos.energies <= xmax + 0.05)
        plot_data = {'mask': mask, 'energies': dos.energies}
        spins = dos.densities.keys()
        ymax = 0

        if plot_total:
            lines = []
            tdos = {'label': 'Total DOS', 'dens': dos.densities,
                    'colour': 'k', 'alpha': 0.15}
            # subplot data formatted as a list of lists of dicts, with each list
            # of dicts being plotted on a seperate graph, if only one list then
            # solo plot
            lines.append([tdos])
            dmax = max([max(d[mask]) for d in dos.densities.values()])
            ymax = dmax if dmax > ymax else ymax
        elif not subplot:
            lines = [[]]  # need a blank list to add lines into

        # TODO: Fix broken behaviour if plot_total is off
        cutoff = (legend_cutoff / 100.) * (ymax / 1.05)

        for el, el_pdos in pdos.iteritems():
            el_lines = []
            for orb in sort_orbitals(el_pdos):
                dmax = max([max(d[mask])
                            for d in el_pdos[orb].densities.values()])
                ymax = dmax if dmax > ymax else ymax
                label = None if dmax < cutoff else '{} ({})'.format(el, orb)
                colour = get_colour_for_element_and_orbital(el, orb, colours)
                el_lines.append({'label': label, 'alpha': 0.25,
                                 'colour': colour,
                                 'dens': el_pdos[orb].densities})
            if subplot:
                lines.append(el_lines)
            else:
                lines[0].extend(el_lines)

        ymax = ymax * empty_space / yscale
        ymin = 0 if len(spins) == 1 else -ymax
        plot_data.update({'lines': lines, 'ymax': ymax, 'ymin': ymin})
        return plot_data

    def get_plot(self, subplot=False, width=6.6, height=5, xmin=-6, xmax=6,
                 yscale=1, colours=None, plot_total=True, legend_on=True,
                 num_columns=2, legend_frame_on=False, legend_cutoff=3, dpi=400,
                 plt=None):
        """Plot the density of states either using matplotlib or xmgrace.

        Args:
            subplot (bool): Split the plot up into separate plots for each
                element.
            width (float): The width of the graph.
            height (float): The height of the graph.
            xmin (float): The minimum energy to plot.
            xmax (float): The maximum energy to plot.
            yscale (dict): Scaling factor for the y-axis.
            colours (dict): Specify custom colours as {'Element': colour} where
                colour is a hex number.
            plot_total (bool): Whether or not to plot total DOS.
            legend_on (bool): Whether or not to plot the graph legend.
            num_columns (int): The number of columns in the legend.
            legend_frame_on (bool): Whether or not to plot the graph legend
                frame.
            legend_cutoff (int): The cut-off (in % of maximum DOS plotted) for a
                elemental/orbital DOS label to appear in the legend.
            dpi (int): The dots-per-inch (pixel density) for the image.
            plt (pyplot object): Matplotlib pyplot object to use for plotting.

        Returns:
            matplotlib pyplot object.
        """
        plot_data = self.dos_plot_data(yscale=yscale, xmin=xmin, xmax=xmax,
                                       colours=colours, plot_total=plot_total,
                                       legend_cutoff=legend_cutoff,
                                       subplot=subplot)

        if subplot:
            nplots = len(plot_data['lines']) + 1
            plt = pretty_subplot(nplots, width=width, height=height, dpi=dpi,
                                 plt=plt)
        else:
            plt = pretty_plot(width=width, height=height, dpi=dpi, plt=plt)

        mask = plot_data['mask']
        energies = plot_data['energies'][mask]
        fig = plt.gcf()
        lines = plot_data['lines']
        spins = [Spin.up] if len(lines[0][0]['dens']) == 1 else \
            [Spin.up, Spin.down]
        for i, line_set in enumerate(plot_data['lines']):
            if subplot:
                ax = fig.axes[i]
            else:
                ax = plt.gca()

            for line in line_set:
                for spin in spins:
                    if spin == Spin.up:
                        label = line['label']
                        densities = line['dens'][spin][mask]
                    elif spin == Spin.down:
                        label = ""
                        densities = -line['dens'][spin][mask]
                    ax.fill_between(energies, densities, lw=0,
                                    facecolor=line['colour'],
                                    alpha=line['alpha'])
                    ax.plot(energies, densities, label=label,
                            color=line['colour'], lw=line_width)

            ax.set_ylim(plot_data['ymin'], plot_data['ymax'])
            ax.set_xlim(xmin, xmax)

            ax.tick_params(axis='x', which='both', top='off')
            ax.tick_params(axis='y', which='both', labelleft='off',
                           labelright='off', left='off', right='off')

            loc = 'upper right' if subplot else 'best'
            ncol = 1 if subplot else num_columns
            # TODO: set size of line in legend (as current mpl overwrite old
            # default)
            if legend_on:
                ax.legend(loc=loc, frameon=legend_frame_on, ncol=ncol,
                          prop={'size': label_size - 3})

        # no add axis labels and sort out ticks
        if subplot:
            fig.text(0.08, 0.5, 'Arb.units', fontsize=label_size, ha='center',
                     va='center', rotation='vertical')
            ax.set_xlabel('Energy (eV)', fontsize=label_size)
            fig.subplots_adjust(hspace=0)
            plt.setp([a.get_xticklabels() for a in fig.axes[:-1]],
                     visible=False)
        else:
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('Arb.units')

        plt.tight_layout()
        return plt


class VBSPlotter(BSPlotter):

    def __init__(self, bs):
        BSPlotter.__init__(self, bs)

    def _maketicks(self, plt):
        """Utility method to add tick marks to a band structure."""
        # set y-ticks
        ax = plt.gca()
        majorLocator = MultipleLocator(2)
        minorLocator = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
        plt.tick_params(axis='both', which='major', labelsize=19)
        plt.tick_params(axis='both', which='minor', labelsize=19)

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
        ax.xaxis.grid(True, c='k', ls='-', lw=vaspy.misc.plotting._linewidth)

    def get_plot(self, zero_to_efermi=True, ymin=-6, ymax=6,
                 width=5, height=5, vbm_cbm_marker=False, plt=None,
                 dos_data=None):
        """
        Get a matplotlib object for the bandstructure plot.
        Blue lines are up spin, red lines are down
        spin.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
        """
        if not plt:
            plt = vaspy.misc.plotting.pretty_plot(width, height, plt=plt)

        band_linewidth = 2

        data = self.bs_plot_data(zero_to_efermi)
        dists = data['distances']
        eners = data['energy']
        nkpts = len(dists[0])
        # nd is branch index, nb is branch index, nk is kpoint index
        for nd, nb in itertools.product(range(len(data['distances'])),
                                        range(self._nb_bands)):
            # colour valence bands blue and conduction bands orange in semiconds
            if (self._bs.is_spin_polarized or self._bs.is_metal() or
                    nb <= self._bs.get_vbm()['band_index'][Spin.up][0] + 1):
                c = '#3953A4'
            else:
                c = '#FAA316'

            e = [eners[nd][str(Spin.up)][nb][nk] for nk in range(nkpts)]
            plt.plot(dists[nd], e, ls='-', c=c, linewidth=band_linewidth)
            if self._bs.is_spin_polarized:
                e = [eners[nd][str(Spin.down)][nb][nk] for nk in range(nkpts)]
                plt.plot(dists[nd], e, 'r--', linewidth=band_linewidth)

        self._maketicks(plt)

        plt.ylabel('Energy (eV)')

        # draw line at Fermi level if not zeroing to e-Fermi
        if not zero_to_efermi:
            ef = self._bs.efermi
            plt.axhline(ef, linewidth=2, color='k')

        # set x and y limits
        plt.xlim(0, data['distances'][-1][-1])
        if self._bs.is_metal() and not zero_to_efermi:
            plt.ylim(self._bs.efermi + ymin, self._bs.efermi + ymax)
        else:
            plt.ylim(ymin, ymax)

        if vbm_cbm_marker:
            for cbm in data['cbm']:
                plt.scatter(cbm[0], cbm[1], color='#D93B2B', marker='o', s=100)
            for vbm in data['vbm']:
                plt.scatter(vbm[0], vbm[1], color='#0DB14B', marker='o', s=100)

        # keep correct aspect ratio square
        ax = plt.gca()
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))

        plt.tight_layout()
        return plt


def get_colour_for_element_and_orbital(element, orbital, colours=None):
    """Select a colour for a particular elemental orbital.

    If that element is not specified in the colours dictionary, a random colour
    will be generated based on the list of 22 colours of maximum contast:
    http://www.iscc.org/pdf/PC54_1724_001.pdf

    Args:
        element (str): The element to select a colour for.
        orbital (str): The orbital.
        colours (dict): A dict of {'Element': {'orb': colour}}, where colour
            is a hex number.

    Returns:
        A colour as either the colour name, hex code, or list of 3 floats.
    """
    try:
        return colours.get(element, orbital)
    except (configparser.NoSectionError, configparser.NoOptionError):
        return col_cycle.next()
