ScriptsForVASP
==============

Making life easier using scripting languages such as bash and python to facilitate multiple job preparation, submission and analysis.

For single-element runs
---------------------
* Prepare.sh entest / lctest / c11 cubic / ... - to create a list of runs
* Fire.sh entest / lctest / c11 / ... - to submit a list of jobs
* Display.sh entest / lctest / c11 / ... - to do the analysis of a list of runs and determine physical quantities

* Prep-fire.sh - a wrapper on Prepare.sh and Fire.sh
* Prep-fast.sh - a fast way to create a folder with necessary input files

* Elastic.sh prep-fire / disp-solve / solve - to do a set of elastic constant determination runs, built on Prepare.sh, Fire.sh and Display.sh
* Dos-bs.sh scrun / dosrun / bsrun / plot-tdos / plot-ldos / plot-bs - to carry out self-consistent runs, density of states runs, band structure runs, to plot the total and local density of states, and to plot the band structure
* Bader.sh prerun /bader - to perform Bader charge transfer analysis

For multi-element runs
---------------------
* ElementRun.sh - an interactive script to execute the above routines for multiple elements

Useful scripts
---------------------
* Cellinfo.sh - to show the important info from OUTCAR
* Swap-file-name.sh - to exchange the names of two files
* Delfiles.sh - to delete the output files and start over for failed runs
* Re-relax.sh - to copy CONTCAR to POSCAR and start over relaxation

General routines
---------------------

1. Create a directory for a specific compound, say PtN. This is the top working directory
2. Go into the directory and create folder INPUT, putting INCAR, KPOINTS, POSCAR, concatenated POTCAR and qsub.parallel into it
3. Customize the input files to prepare the run, including reserving the strings that need to be replaced as @R@, @N@, etc.
4. Go back to the top working directory. Use Prepare.sh xxx yyy to perform a set of runs, or Prep-fast.sh xxx to simply create a folder xxx with proper substitutions
5. Use Fire.sh xxx to execute a set of qsub batch commands, or simply go into the folder prepared by Prep-fast.sh and qsub qsub.parallel
6. Using Display.sh xxx [M/P (lctest only)] to get the (step, energy) pairs written into xxx_output.txt and other instructive information. For lctest and elastic constant runs the fitting results and plot figures are also written
7. Rename the folder created during the process under the top working directory to whatever looks neat and tidy, i.e. labeling them with numbers

To perform a complete set of elastic constant runs
---------------------

1. Under the top working directory, execute Elastic.sh xxx prep-fire, i.e. cubic
2. Elastic.sh xxx disp-solve after every batched job is finished
3. Go check the elastic/elastic_output.txt for the results

To perform runs for multiple elements
---------------------

1. Create a directory, i.e. TMN
2. Create a folder INPUT, putting in necessary input files with necessary @R@ and @N@
2. Go into the directory and execute Element-run.sh. Do as suggested by the script
