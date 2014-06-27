Scripts structure
=================

Python scripts
--------------

These can also be used alone, and **can be readily imported and played with in a Python interpreter**.

* ``Plot_tdos.py`` - to plot total density of states, spin states are supported (above or below x axis)
* ``Plot_ldos.py`` - to plot local projected density of states of two ions (above or below x axis), spin states are supported (in different files)
* ``Plot_bs.py`` - to plot the band structure along certain chosen points in the 1st BZ
* ``Plot_cohp.py`` - to plot the pCOHP from LOBSTER runs
* ``Strain_applier.py`` - to add strain to the cell lattice vectors in POSCAR for a series of runs, determining the value of a combination of elastic constants
* ``Ions_rotator.py`` - to rotate ions in selected lines in POSCAR at a certain point, around a certain axis, to a certain angle

Bash convenient scripts
-----------------------

* ``Cell_info.sh [rwigs]`` - to show the important info from OUTCAR (rwigs is optional, and displays the RWIGS values for different species to 100 % fill the cell)
* ``Del_files.sh [XXX]`` - to delete the output files for failed runs (XXX is optional if one stays in the job directory)
* ``Swap_file_names.sh AAA BBB`` - to exchange the names of two files AAA and BBB

Bash Job submission scripts
---------------------------

* ``Prepare.sh XXX`` - a fast way to create a directory XXX with necessary input files
* ``Sequence_test.sh TEST_TYPE [OPTIONS]`` - to create a series of runs, performing encut test | kpoints test | lattice constant test | elastic constant test | ...
* ``Sequence_disp.sh TEST_TYPE [OPTIONS]`` - to do the post-analysis of a series of runs described above and determine physical quantities
* ``Elastic.sh test | disp-solve | solve`` - to do a full set of independent elastic constant determination runs
* ``Electronic.sh scrun | dosrun | bsrun | lobster kp | lobster test | lobster analysis``
  - to carry out self-consistent run, density of states run, band structure run, LOBSTER pCOHP run
* ``Bader.sh test | analysis`` - to perform Bader charge transfer analysis

Bash multiple element species
-----------------------------

* ``Element_run.sh`` - an interactive script to execute the above routines for multiple species

Supportive scripts
------------------

**These are used by Bash scripts and better left unvisited.**

* ``_qsub_replacer.sh`` - to do some replacement to the job submission file so that it corresponds to the right working directory. Machine dependent
* ``_sequence_fit_plot.py`` - to perform different fitting techniques, including Birch-Murnaghan equation of state and polynomial
* ``_elastic_solver.py`` - to solve the independent elastic constants from the same number of combination runs with linear algebra
* ``_cell_info_solver.py`` - to obtain the RWIGS values for different species to 100 % fill the cell

**others** contains some interesting physics or computational experiments.