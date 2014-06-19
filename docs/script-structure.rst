Scripts structure
=================

Basic level
-----------

* ``Fast-prep.sh XXX`` - a fast way to create a directory XXX with necessary input files
* ``Prepare.sh TEST_TYPE [OPTIONS]`` - to create a series of runs, performing encut test | kpoints test | lattice constant test | elastic constant test | ...
* ``Fire.sh TEST_TYPE [OPTIONS]`` - to submit a series of jobs described above
* ``Display.sh TEST_TYPE [OPTIONS]`` - to do the post-analysis of a series of runs described above and determine physical quantities
* ``Prep-fire.sh TEST_TYPE [OPTIONS]`` - a wrapper on Prepare.sh and Fire.sh

Higher level
------------

* ``Elastic.sh prep-fire | disp-solve | solve`` - to do a full set of independent elastic constant determination runs.
* ``Electronic.sh scrun | dosrun | bsrun | lobster-kp | lobster-prerun | lobster | plot-tdos | plot-ldos | plot-bs | plot-cohp``
  - to carry out self-consistent runs, density of states runs, band structure runs, LOBSTER pCOHP runs
  and to plot the total and local density of states, band structure, and pCOHP
* ``Bader.sh prerun | bader`` - to perform Bader charge transfer analysis

For multiple element species
----------------------------

* ``ElementRun.sh`` - an interactive script to execute the above routines for multiple species

For convenience
---------------

* ``Cellinfo.sh [rwigs]`` - to show the important info from OUTCAR (rwigs is optional, and displays the RWIGS values for different species to 100 % fill the cell)
* ``Delfiles.sh [XXX]`` - to delete the output files for failed runs (XXX is optional if one stays in the job directory)
* ``Re-relax.sh [XXX]`` - to copy CONTCAR to POSCAR and start over relaxation (XXX is optional if one stays in the job directory)
* ``Swap-files-name.sh AAA BBB`` - to exchange the names of two files AAA and BBB

Python scripts
--------------

These can also be used alone, and **can be readily imported and played with**.

* ``Plot_tdos.py`` - to plot total density of states, spin states are supported (above or below x axis)
* ``Plot_ldos.py`` - to plot local projected density of states of two ions (above or below x axis), spin states are supported (in different files)
* ``Plot_bs.py`` - to plot the band structure along certain chosen points in the 1st BZ
* ``Plot_cohp.py`` - to plot the pCOHP from LOBSTER runs
* ``Ions-rotator.py`` - to rotate ions in selected lines in POSCAR at a certain point, around a certain axis, to a certain angle

These are better left unvisited.

* ``_prepare_strain.py`` - to add strain to the cell lattice vectors in POSCAR for a series of runs, determining the value of a combination of elastic constants
* ``_display_fit.py`` - to perform different fitting techniques, including Birch-Murnaghan equation of state and polynomial
* ``_elastic_solver.py`` - to solve the independent elastic constants from the same number of combination runs with linear algebra
* ``_cellinfo_solver.py`` - to obtain the RWIGS values for different species to 100 % fill the cell