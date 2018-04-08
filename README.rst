Vaspy
=====

Vaspy is a Python toolkit for plotting and analysing VASP calculations. The main features include:

 1. Extensive framework for generating high-symmetry k-point paths.
 2. Plotting scripts for electronic and phonon band structures, density of states, and optical absorption.
 3. Analysis scripts to calculate band edge effective masses.

Vaspy is built on pymatgen, meaning there is no effort required to read and write typical VASP files.
Vaspy is free to use, however, we welcome your help in improving and extending the
package with your own contributions.

Warning: There are probably still some bugs. If you think you've find any,
please report them on the Issue Tracker.

Usage
-----

Primarily, vaspy is intended to be used via the command-line, however, a fully-documented
python API is also provided. The package documentation is still a WIP, however the built-in
help (`` -h``) option for each command provides a summary of the available options.

Currently the scripts provided by vaspy are:
 - ``vaspy-kgen``: For generating VASP KPOINTS files along high-symmetry k-point paths.
 - ``vaspy-bandplot``: For plotting publication-ready electronic band structure diagrams.
 - ``vaspy-dosplot``: For plotting publication-ready electronic density of state diagrams.
 - ``vaspy-optplot``: For plotting publication-ready optical absorption plots.
 - ``vaspy-phonon-bandplot``: For plotting publication-ready phonon band structure diagrams.
 - ``vaspy-bandstats``: For calculating electron and hole effective masses from a band structure.

Installation
------------

We recommend installation from source with Pip, this will automatically install any dependencies::

  pip3 install --user .

Developer installation
~~~~~~~~~~~~~~~~~~~~~~

Developers may prefer to install using ``pip3 install --user -e .`` which
creates an "editable" local installation. Instead of copying files,
this creates links to the source folder so that that tweaks to the
code in your source folder will be immediately reflected on the PATH.

Tests
~~~~~

To ensure the code has been installed correctly, the unittests can be run using::

  python -m unittest discover tests

Requirements
------------

Vaspy is currently compatible with Python 3.4+ and relies on a number of
open-source python packages, specifically:

- Pymatgen
- Numpy
- Scipy
- Matplotlib
- Spglib
- Phonopy
- H5py (optional dependency for phonon plotting features)

Vaspy uses Pip and setuptools for installation. You *probably* already
have this; if not, your GNU/Linux package manager will be able to oblige
with a package named something like ``python-setuptools``. On Max OSX
the Python distributed with `Homebrew <http://brew.sh>`_. includes
setuptools and Pip.

License
-------

Vaspy is made available under the MIT License.
