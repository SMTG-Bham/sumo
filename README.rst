Sumo
====

.. image:: https://travis-ci.org/SMTG-UCL/sumo.svg?branch=master
    :target: https://travis-ci.org/SMTG-UCL/sumo

.. image:: https://readthedocs.org/projects/sumo/badge/?version=latest
    :target: http://sumo.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Sumo is a Python toolkit for plotting and analysis of ab initio
calculation data. The main features include:

1. An extensive framework for generating high-symmetry k-point paths.
2. Plotting scripts for electronic and phonon band structures, density
   of states, and optical absorption diagrams.
3. Analysis scripts to calculate parabolic and non-parabolic band
   effective masses.

The code currently only supports VASP calculations, however, we plan to
add support for additional codes in future releases.

Sumo is free to use, however, we ask that you cite the code if you use
it in your research.

Warning: There are probably still some bugs. If you think you've found
one, please report it on the `Issue Tracker
<https://github.com/SMTG-UCL/sumo/issues>`_. We welcome your help in
improving and extending the package with your own contributions.


Usage
-----

Sumo is intended to be used via the command-line, however, a
fully-documented python API is also provided. A manual, including
tutorials and API documentation, is `available online
<http://sumo.readthedocs.io/en/latest/>`_. Additionally, the built-in
help (``-h``) option for each command provides a summary of the
available options.

Currently, the scripts provided by sumo are:

- ``sumo-kgen``: For generating VASP KPOINTS files along high-symmetry
  k-point paths.
- ``sumo-bandplot``: For plotting publication-ready electronic band
  structure diagrams.
- ``sumo-dosplot``: For plotting publication-ready electronic density of
  states diagrams.
- ``sumo-optplot``: For plotting publication-ready optical absorption
  diagrams.
- ``sumo-phonon-bandplot``: For plotting publication-ready phonon band
  structure diagrams.
- ``sumo-bandstats``: For calculating electron and hole effective masses
  from a band structure.

A guide to using each command can be found on the
`Tutorial page <http://sumo.readthedocs.io/en/latest/tutorials.html>`_.

For a preview of the functionality of sumo, see the
`Gallery <http://sumo.readthedocs.io/en/latest/gallery.html>`_.


Installation
------------

We recommend installation from source with Pip, this will automatically
install any dependencies:

.. code-block:: bash

    pip3 install --user sumo

To build the documentation, download the package source and install with
extra dependecies:

.. code-block:: bash

    pip3 install --user .[docs]
    cd docs
    make html


Developer installation
~~~~~~~~~~~~~~~~~~~~~~

Developers may prefer to install using ``pip3 install --user -e .``
which creates an "editable" local installation. Instead of copying files,
this creates links to the source folder so that that tweaks to the
code in your source folder will be immediately reflected on the PATH.


Tests
~~~~~

To ensure the code has been installed correctly, the unittests can be
run (from the root directory of the project) using::

  python -m unittest discover tests


Requirements
------------

Sumo is currently compatible with Python 3.5+ and relies on a number of
open-source python packages, specifically:

- `Pymatgen <http://pymatgen.org>`_
- `Numpy <http://www.numpy.org>`_
- `Scipy <https://www.scipy.org>`_
- `Matplotlib <https://matplotlib.org>`_
- `Spglib <https://atztogo.github.io/spglib/>`_
- `Phonopy <https://atztogo.github.io/phonopy/>`_
- `SeeK-path <https://github.com/giovannipizzi/seekpath>`_
- `H5py <https://www.h5py.org>`_

Sumo uses Pip and setuptools for installation. You *probably* already
have this; if not, your GNU/Linux package manager will be able to oblige
with a package named something like ``python-setuptools``. On Max OSX
the Python distributed with `Homebrew <http://brew.sh>`_. includes
setuptools and Pip.


License
-------

Sumo is made available under the MIT License.
