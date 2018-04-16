Sumo
====

Sumo is a Python toolkit for plotting and analysing VASP calculations. The main features include:

1. An extensive framework for generating high-symmetry k-point paths.
2. Plotting scripts for electronic and phonon band structures, density of states, and optical absorption.
3. Analysis scripts to calculate band edge effective masses.

Sumo is built on pymatgen, meaning there is no effort required to read and write typical VASP files.
Sumo is free to use, however, we welcome your help in improving and extending the
package with your own contributions.

Warning: There are probably still some bugs. If you think you've find any,
please report them on the Issue Tracker.

Usage
-----

Primarily, sumo is intended to be used via the command-line, however, a fully-documented
python API is also provided.
A full manual, including tutorials and API documentation,
is available online at `readthedocs.io <http://sumo.readthedocs.io/en/latest/>`__.
Additionally, the built-in help (``-h``) option for each command provides a
summary of the available options.

Currently the scripts provided by sumo are:

- ``sumo-kgen``: For generating VASP KPOINTS files along high-symmetry k-point paths.
- ``sumo-bandplot``: For plotting publication-ready electronic band structure diagrams.
- ``sumo-dosplot``: For plotting publication-ready electronic density of state diagrams.
- ``sumo-optplot``: For plotting publication-ready optical absorption plots.
- ``sumo-phonon-bandplot``: For plotting publication-ready phonon band structure diagrams.
- ``sumo-bandstats``: For calculating electron and hole effective masses from a band structure.

Installation
------------

We recommend installation from source with Pip, this will automatically install any dependencies:

.. code-block:: bash

    pip3 install --user .

To build the documentation, install the package with extra dependecies

.. code-block:: bash

    pip3 install --user .[docs]
    cd docs
    make html

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

Sumo is currently compatible with Python 3.4+ and relies on a number of
open-source python packages, specifically:

- Pymatgen
- Numpy
- Scipy
- Matplotlib
- Spglib
- Phonopy
- H5py (optional dependency for phonon plotting features)

Sumo uses Pip and setuptools for installation. You *probably* already
have this; if not, your GNU/Linux package manager will be able to oblige
with a package named something like ``python-setuptools``. On Max OSX
the Python distributed with `Homebrew <http://brew.sh>`_. includes
setuptools and Pip.

License
-------

Sumo is made available under the MIT License.
