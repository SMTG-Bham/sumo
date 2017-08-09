from matplotlib.ticker import MultipleLocator
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

# TODO:
#  - set font to whitney


def plot_band_structure(band_structure, ymin=-6, ymax=6, height=6, width=6,
                        plt_format='mpl', vbm_cbm_marker=False):
    plotter = BSPlotter(band_structure)
    plt = plotter.get_plot(vbm_cbm_marker=vbm_cbm_marker)
    fig = plt.gcf()
    ax = plt.gca()

    # set axis properties
    majorLocator = MultipleLocator(2)
    minorLocator = MultipleLocator(1)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    plt.tick_params(axis='both', which='major', labelsize=19)
    plt.tick_params(axis='both', which='minor', labelsize=19)
    plt.ylabel("Energy (eV)", size=20)
    plt.xlabel("")
    plt.ylim(ymin, ymax)

    # re-colour vbm and cbm
    vbm = band_structure.get_vbm()["band_index"][Spin.up][0]
    lines = ax.get_lines()
    for i in range(band_structure.nb_bands):
        c = '#3953A4' if i <= vbm + 1 else '#FAA316'
        plt.setp(lines[i], color=c, linewidth=1.5)

    if vbm_cbm_marker:
        # re-colour the vbm and cbm indicators
        ax.get_children()[1].set_color('#D93B2B')
        ax.get_children()[0].set_color('#0DB14B')

    # make the figure and plot are square
    fig.set_size_inches(width, height)
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_aspect((x1-x0)/(y1-y0))
    return plt
