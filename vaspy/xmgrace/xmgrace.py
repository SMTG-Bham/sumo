import os, sys, shutil

from vaspy.vasp_input import *
from vaspy.vasp_output import *
from vaspy.math_func import *

class XmgraceDocument:

    """
    Generic class to quickly plot band structures/dos

    Args:
        band: (bool) of whether to plot as band or doscar
    """
    def __init__(self, band, height=792, width=612):
        self._band = band
        self._graphs = []
        self._height = height
        self._width = width
        self._graph_sizes = []
        self._total_kpoints = 0



    def add_graph(self, graph):
        self._graphs.append(graph)
        self._total_kpoints += graph.num_kpoints

        if len(self._graphs) == 1:
            self._graph_sizes.append([0.15, 0.15, 0.85, 0.85])
        else:
            self._graph_sizes = []
            x_cum = 0.15
            avail_space = 0.7 - (len(self._graphs) - 1) * 0.082059
            for graph in self._graphs:
                perc_width = float(graph.num_kpoints)/float(self._total_kpoints)
                width = perc_width * avail_space
                self._graph_sizes.append([x_cum, 0.15, x_cum+width, 0.85])
                x_cum += width
                x_cum += 0.082059


    def to_file(self, filename):
        if self._band:
            base = os.path.dirname(os.path.realpath(__file__)) + '/misc/base.agr'
        shutil.copy(base, filename)

        with open(filename, "r+") as f:
            content = f.read()
            f.seek(0, 0)
            f.write("@page size %d, %d\n%s" % (self._height, self._width, content))

            graph_id = 0
            for graph, size in zip(self._graphs, self._graph_sizes):
                f.write(graph.str_graph_settings(graph_id, size))
                graph_id += 1

            graph_id = 0
            for graph in self._graphs:
                f.write(graph.str_data(graph_id))
                graph_id += 1
                #print graph_id



class XmgraceGraph:

    def __init__(self, pos=0, xy_data=None, num_kpoints=None,
                 data_settings=None, axis_labels=None, emin=-6, emax=6):
        if xy_data is None:
            xy_data = []
        if num_kpoints is None:
            num_kpoints = []
        if data_settings is None:
            data_settings = []
        if axis_labels is None:
            axis_labels = []

        self._data = xy_data
        self.num_kpoints = num_kpoints
        self._data_settings = data_settings
        self._axis_labels = axis_labels
        if pos == 0:
            self._pos = "normal"
        if pos == 1:
            self._pos = "opposite"
        if pos == 2:
            self._pos = "none"

        self.emin = emin
        self.emax = emax

    """
    Args:
        data: a list of [x, y] values
        color: the line color (4 = blue, 11 = yellow)
        line_width: the width of the line
        label: label for key
        comment: a comment to aid identification of the data set
    """
    def add_xy(self, data, line_color, line_width, label="", comment=""):
        self._data.append(data)
        self._data_settings.append({'line_width': line_width, 'label': label,
            'comment': comment, 'line_color': line_color})

    def set_axis_labels(self, labels):
        self._axis_labels = labels

    def str_graph_settings(self, graph_id, size):
        s = "@g%d on\n" % graph_id
        s += "@g%d on\n" % graph_id
        s += "@g%d hidden false\n" % graph_id
        s += "@g%d type XY\n" % graph_id
        s += "@g%d stacked false\n" % graph_id
        s += "@g%d bar hgap 0.000000\n" % graph_id
        s += "@g%d fixedpoint off\n" % graph_id
        s += "@g%d fixedpoint type 0\n" % graph_id
        s += "@g%d fixedpoint xy 0.000000, 0.000000\n" % graph_id
        s += "@g%d fixedpoint format general general\n" % graph_id
        s += "@g%d fixedpoint prec 6, 6\n" % graph_id
        s += "@with g%d\n" % graph_id
        s += "@    world 1, %d, %d, %d\n" % (self.emin, self.num_kpoints, self.emax)
        s += "@    stack world 0, 0, 0, 0\n"
        s += "@    znorm 1\n"
        s += "@    view %.6f, %.6f, %.6f, %.6f\n" % (size[0], size[1], size[2], size[3])
        settings_template = base = os.path.dirname(os.path.realpath(__file__)) + '/misc/graph.agr'
        with open(settings_template) as f:
            lines = f.readlines()

        s += "".join(lines)

        s += self.str_band_labels()

        settings_id = 0
        for settings in self._data_settings:
            s += self.str_data_settings(settings, settings_id)
            settings_id += 1

        return s

    def str_band_labels(self):
        s = "@    xaxis  tick spec type both\n"
        s += "@    xaxis  tick spec %d\n" % len(self._axis_labels)
        s += "@    yaxis  ticklabel place %s\n" % self._pos
        title = "Energy (eV)" if self._pos == "normal" else ""
        s += "@    yaxis  label \"%s\"\n" % title
        label_id = 0
        for label, start in self._axis_labels:
            s += "@    xaxis  tick major %d, %d\n" % (label_id, start)
            s += "@    xaxis  ticklabel %d, \"%s\"\n" % (label_id, label)
            label_id += 1
        return s

    def str_data_settings(self, settings, settings_id):
        settings_template = base = os.path.dirname(os.path.realpath(__file__)) + '/misc/data_settings'
        with open(settings_template) as f:
            lines = f.readlines()

        prepend = "@    s%d " % settings_id
        prepended_lines = [prepend + s for s in lines]
        s = "".join(prepended_lines)
        s += "@    s%d line linewidth %.1f\n" % (settings_id, settings['line_width'])
        s += "@    s%d line color %d\n" % (settings_id, settings['line_color'])
        s += "@    s%d comment \"%s\"\n" % (settings_id, settings['comment'])
        s += "@    s%d legend \"%s\"\n" % (settings_id, settings['label'])

        return s

    def str_data(self, graph_id):
        s = ""
        data_id = 0
        for data in self._data:
            s += "@target G%d.S%d\n" % (graph_id, data_id)
            s += "@type xy\n"
            for x, y in data:
                s += "%s %s\n" % (str(x), str(y))
            s += "&\n"
            data_id += 1

        return s
