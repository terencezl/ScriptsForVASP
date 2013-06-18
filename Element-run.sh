#!/bin/bash
# Use: In the TMN working directory
# Change element.dat first

element_list=$(echo $(cat INPUT_ELEMENT/element.dat))
echo "Tell me what you wish to do: prepare / lctest / display-M / M-to-P / display-P / elastic / solve / alternative / solve-again / scrun / dosrun / plotdos / elastic-scrun / elastic-dosrun / elastic-plotdos / bader-prerun / bader"
read -r test_type
echo "What lattice type are you thinking of? zinc-blende / rocksalt / CsCl /..."
read -r cryst_struct

# create INPUT folder in the selected cryst structure of each element
if [ $test_type == prepare ]; then
    for n in $element_list
    do
        mkdir -p $n/$cryst_struct/INPUT
        cp -r INPUT_ELEMENT/{INCAR,KPOINTS,POSCAR_$cryst_struct,qsub.parallel} $n/$cryst_struct/INPUT
        # select LDA or GGA, concatenate metal and nitrogen
        cat INPUT_ELEMENT/POTCAR_GGA/{POTCAR_$n,POTCAR_N} > $n/$cryst_struct/INPUT/POTCAR
#        cat INPUT_ELEMENT/POTCAR_LDA/{POTCAR_Na,POTCAR_$n} > $n/$cryst_struct/INPUT/POTCAR
        cd $n/$cryst_struct/INPUT
        mv POSCAR_$cryst_struct POSCAR
        sed -i "s/@N@/$n C ($cryst_struct)/g" INCAR
        sed -i "s/@N@/$n C ($cryst_struct)/g" POSCAR
        cd ../../..
    done

# basic lctest in the "top working directory"
elif [ $test_type == lctest ]; then
    echo "Enter the starting scaling factor, ending scaling factor and the step, separated by a space."
    read -r start end step
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Prep-fire.sh lctest $start $end $step
        cd ../..
    done

# fit the lctest curve to Murnaghan eqn of state
elif [ $test_type == display-M ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Display.sh lctest M
        cd ../..
    done

# rename lctest to lctest_M, grep the lattice constant and do a 5-point close local lctest
elif [ $test_type == M-to-P ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        M-to-P.sh
        cd ../..
    done

# fit the lctest curve to the 2nd polynomial
elif [ $test_type == display-P ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Display.sh lctest P
        cd ../..
    done

# prep and fire the elastic constant runs
elif [ $test_type == elastic ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Elastic.sh $cryst_sys prep-fire
        cd ../..
    done

# display the fitting results and solve the elastic constants
elif [ $test_type == solve ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Elastic.sh $cryst_sys disp-solve
        cd ../..
    done

# prep and fire a different strain
elif [ $test_type == alternative ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    if [ $cryst_sys == cubic ]; then
        alternative='2c11+2c12+c44'
    fi
    for n in $element_list
    do
        cd $n/$cryst_struct/elastic || exit 1
        Prep-fire.sh $alternative $cryst_sys
        cd ../../..
    done

# solve the elastic constants with the new set of strains
elif [ $test_type == solve-again ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    if [ $cryst_sys == cubic ]; then
        alternative='2c11+2c12+c44'
    fi
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        (cd elastic; Display.sh $alternative)
        Elastic-solve.sh $cryst_sys alternative
        cd ../..
    done

# scrun - do the self-consistent run to get the CHGCAR
# dosrun - do the dos run in the sc-dos-bs folder, putting the last run into a deeper folder called scrun
elif [ $test_type == scrun ] || [ $test_type == dosrun ] || [ $test_type == plotdos ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Dos-bs.sh $test_type
        cd ../..
    done

# scrun and dos run for different elastic constants
elif [ $test_type == elastic-scrun ] || [ $test_type == elastic-dosrun ] || [ $test_type == elastic-plotdos ]; then
    test_type=${test_type#elastic-}
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Elastic-dos.sh $cryst_sys $test_type
        cd ../..
    done
    
elif [ $test_type == bader-prerun ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Prep-fast.sh bader
        scaling_factor=$(grep "Equilibrium scaling factor is" lctest/lctest_output.txt | head -1 | awk '{print $5}')
        cd bader
        sed -i "s/@R@/$scaling_factor/g" POSCAR
        echo 'LAECHG = .TRUE. ' >> INCAR
        sed -i 's/LCHARG =/#LCHARG =/g' INCAR
        qsub qsub.parallel
        cd ../../..
    done
    
elif [ $test_type == bader ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct/bader || exit 1
        chgsum.pl AECCAR0 AECCAR2
        bader CHGCAR -ref CHGCAR_sum
        cd ../../..
    done
    grep '' */*/*/ACF*

# template
elif [ $test_type == copy ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        cp sc-dos-bs/LDOS* ../..
        cd ../..
    done
    
fi
