Main Task Scripts:

::::::Multi-element runs::::::
ElementRun.sh - Interactive

::::::Elastic set runs::::::
Elastic-dos.sh cubic/... prerun / scrun / dosrun / plotdos
Elastic.sh cubic/... prep-fire / disp-solve
Elastic-solve.sh cubic/... O / P
Elastic-moduli.sh cubic/...
    _Elastic-solve_solver.py cubic/... volumn_of_primitive_cell input_data(in the form of a python list) original / alternative

::::::Single batch runs::::::
Dos-bs.sh scrun / dosrun / plotdos
MtoP.sh
Display.sh entest/lctest/c11... M / P (for lctest)
    _Display-fit.py entest / lctest / c11... line_from(python style) line_to(python style) M / P (for lctest)
Prep-fire.sh entest / lctest / c11/... [cryst_sys OR from to step (lctest only)]
Prepare.sh entest / lctest / c11/... cubic/...
    _Prepare-strain.py c11/... cubic/... delta
Fire.sh entest / lctest / c11/... [from to step (lctest only)]

----------------------------------------------------

Auxiliary Scripts:

Cellinfo.sh [rwigs] - Go into the deepest VASP working directory and execute
    _Cellinfo-solver.py task(rwigs) volumn_of_primitive_cell input_data(in the form of a python list)

----------------------------------------------------

Small Quick Scripts:

Delfiles.sh test_type(existent subfolder of the current directory)
PrepRerun.sh test_type(existent subfolder of the current directory) - Clean up some output files that may influence the next run 

####################################################################################

General Uses:

1. Create a directory for a specific compound, say PtN. This is the top working directory
2. Go into the directory and create folder INPUT, putting INCAR, KPOINTS, POSCAR, concatenated POTCAR and qsub.parallel into it
3. Customize the input files to prepare the run, including reserving the strings that need to be replaced as @R@, @N@, etc.
4. Go back to the top working directory. Use Prepare.sh xxx yyy to perform a set of runs, or Prep-fast.sh xxx to simply create a folder xxx with proper substitutions
5. Use Fire.sh xxx to execute a set of qsub batch commands, or simply go into the folder prepared by Prep-fast.sh and qsub qsub.parallel
6. Using Display.sh xxx [M/P (lctest only)] to get the (step, energy) pairs written into xxx_output.txt and other instructive information. For lctest and elastic constant runs the fitting results and plot figures are also written
7. Rename the folder created during the process under the top working directory to whatever looks neat and tidy, i.e. labeling them with numbers

----------------------------------------------------

To perform a complete set of elastic constant runs:

1. Under the top working directory, execute Elastic.sh xxx prep-fire, i.e. cubic
2. Elastic.sh xxx disp-solve after every batched job is finished
3. Go check the elastic/elastic_output.txt for the results

----------------------------------------------------

To use a slightly different set of elastic constant equations (i.e. changing only one equation):

1. Perform an original complete run with Elastic.sh xxx prep-fire, Elastic.sh xxx disp-solve
2. Go into the elastic folder to perform an additional run (with a specific combination of elastic consts) with Prepare.sh, Fire.sh and Display.sh
3. Go out of the folder and execute Elastic-solve.sh xxx alternative, the results will show

----------------------------------------------------

To perform a large scale run through different elements:

1. Create a directory, i.e. TMN
2. Create a folder INPUT, putting in necessary input files with necessary @R@ and @N@
2. Go into the directory and execute Element-run.sh. Do as suggested by the script
