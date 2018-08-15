Change Log
==========

[Unreleased]
------------

v1.1.1
-------

Fix bug when installing from Pypi.


v1.1.0
-------

Use matplotlib style sheets for styling plots (@ajjackson & @utf).
Enables plots to be customised based on user settings.

Various bug fixes:

- Fix bug when normalising DOS to Fermi level.
- Fix codacy style issues.
- Plotting style standardised across all plots.

v1.0.10
-------

Add option to align DOS to Fermi level (@shyamd)

Various bug fixes:

- Fix many typos.
- Updates to paper and documentation.

v1.0.9
------

``phonon-bandplot`` now supports combined DOS & band structure plots (Adam Jackson, Arthur Yaud).

Various bug fixes:

- Fix P centered trigonal k-point path.
- Fix ``--symprec`` behaviour in phonon-bandplot.
- Fix orbital projected band structures with branches (Adam Jackson).
- Fix reading Eg from spin-pol calculations (Adam Jackson).

v1.0.8
------

Enhancements by Adam Jackson:

- Add y-label and dos label options for DOS & band plots.
- Cache DOS colours for consistent plots.

Various bug fixes:

- Fixed gaussian broadening of DOS.
- Fixed ``--spg`` option in kgen and phonon-bandplot.
- Fixed default arguments for band structure + dos plotting.
- Added A centered orthorhombic lattice to ``BradCrackKpath``.

v1.0.7
------

Various bug fixes:

- Fixed density option in kgen.
- Fixed phonon-bandplot plotting limits.

v1.0.6
------

Move package data files.

v1.0.5
------

Minor bug fixes.

v1.0.4
------

Minor changes to Pypi config.

v1.0.0
------

Added
~~~~~

- Script files:

  - sumo-kgen
  - sumo-dosplot
  - sumo-bandplot
  - sumo-bandstats
  - sumo-optplot
  - sumo-phonon-bandplot

`[Unreleased] <https://github.com/smtg-ucl/sumo/compare/v1.0.9...HEAD>`_.
`[1.0.8] <https://github.com/smtg-ucl/sumo/compare/v1.0.4...v1.0.9>`_.
