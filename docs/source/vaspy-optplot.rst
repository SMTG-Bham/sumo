vaspy-optplot
==============

``vaspy-optplot`` is a program for generating publication-ready optical absorption
spectra diagrams from VASP calculations. The script supports plotting multiple
spectra simultaneously.

.. contents:: Table of Contents
   :local:
   :backlinks: None

Usage
-----

The full range of options supported by ``vaspy-optplot`` are detailed in the `Command-Line Interface`_ section,
and be can be accessed using the command::

    vaspy-optplot -h

To plot an absorption spectra, simply run the following command in a folder containing a ``vasprun.xml`` or
``vasprun.xml.gz`` file, which has been calculated using ``LOPTICS = .TRUE.``::

    vaspy-optplot

The plot will be written to a file named ``absorption.pdf``, with the raw data written to ``absorption.dat``.

For example, if we run the command in the ``vaspy/tests/data/Cs2SnI6/optics`` directory, the absorption
spectra should look like:

.. image:: figures/absorption_basic.png
   :height: 400px
   :align: center


Basic Options
~~~~~~~~~~~~~

``vaspy-optplot`` automatically searches for a ``vasprun.xml`` file in the current directory.
To specify a particular ``vasprun.xml`` to plot, the ``--filenames`` option can be used.

The height, and width of the graphic, along with the y-axis limits, can be controlled via the
``--width``, ``--height``, ``--ymax``, and ``--ymin`` options.

Additional gaussian broading can be applied using the ``--gaussian`` option. The setting expects a floating
point number as the argument and controls the standard deviation of the broadening applied.


Anisotropic Absorption
~~~~~~~~~~~~~~~~~~~~~~

By default, ``vaspy-optplot`` plots the average optical absorption. The anisotropic contributions
from the x, y and z cartesian directions can be plotted individually using the ``--anisotropic``
option.

For example, if we run the following command in the ``vaspy/tests/data/Cs2SnI6/optics`` directory,
the anisotropic absorption spectra should look like::

    vaspy-optplot --anisotropic

.. image:: figures/absorption_anisotropic.png
   :height: 400px
   :align: center


Plotting Multiple Spectra
~~~~~~~~~~~~~~~~~~~~~~~~~



Things to cover:
- plotting multiple spectra
- band gaps
- gaussian broadening
- ansiotropic plotting

Command-Line Interface
----------------------

.. argparse::
   :module: vaspy.cli.optplot
   :func: _get_parser
   :prog: vaspy-optplot
