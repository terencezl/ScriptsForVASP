#!/bin/bash
# Use: In the TMN working directory
# Change element.dat first

element_list=$(echo $(cat INPUT_ELEMENT/element.dat))
echo "Tell me what you wish to do: prepare / lctest / display_M / lctest_M_to_P / display_P / elastic / solve / alternative / solve_again"
read -r test_type
echo "What lattice type are you thinking of? zinc-blende / rocksalt / CsCl /..."
read -r cryst_struct

# create INPUT folder in the selected cryst structure of each element
if [ $test_type == prepare ]; then
    for n in $element_list
    do
        mkdir -p $n/$cryst_struct/INPUT
        cp -r INPUT_ELEMENT/{INCAR,KPOINTS,POSCAR_$cryst_struct,qsub.parallel} $n/$cryst_struct/INPUT
        cat INPUT_ELEMENT/POTCAR_GGA/{POTCAR_$n,POTCAR_N} > $n/$cryst_struct/INPUT/POTCAR
        cd $n/$cryst_struct/INPUT
        mv POSCAR_$cryst_struct POSCAR
        sed -i "s/@N@/$n N ($cryst_struct)/g" INCAR
        sed -i "s/@N@/$n N ($cryst_struct)/g" POSCAR
        cd ../../..
    done

# basic lctest in the "top working directory"
elif [ $test_type == lctest ]; then
    echo "Enter the starting scaling factor, ending scaling factor and the step, separated by a space."
    read -r start end step
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Prepare.sh lctest $start $end $step
        Fire.sh lctest $start $end $step
        cd ../..
    done

# fit the lctest curve to Murnaghan eqn of state
elif [ $test_type == display_M ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Display.sh lctest M
        cd ../..
    done

# rename lctest to lctest_M, grep the lattice constant and do a 5-point close local lctest
elif [ $test_type == lctest_M_to_P ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        rename lctest lctest_M lctest
        scaling_factor_from_M_fit=$(printf %.2f $(grep "Equilibrium scaling factor is" lctest_M/lctest_output.txt | head -1 | awk '{print $5}'))
        scaling_factor_start=$(echo "$scaling_factor_from_M_fit - 0.03" | bc)
        scaling_factor_end=$(echo "$scaling_factor_from_M_fit + 0.03" | bc)
#        Prepare.sh lctest $scaling_factor_start $scaling_factor_start 0.01
#        Prepare.sh lctest $scaling_factor_end $scaling_factor_end 0.01
#        Fire.sh lctest $scaling_factor_start $scaling_factor_start 0.01
#        Fire.sh lctest $scaling_factor_end $scaling_factor_end 0.01
        Prepare.sh lctest $scaling_factor_start $scaling_factor_end 0.01
        Fire.sh lctest $scaling_factor_start $scaling_factor_end 0.01
        cd ../..
    done

# fit the lctest curve to the 2nd polynomial
elif [ $test_type == display_P ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Display.sh lctest P
        cd ../..
    done

# prep and fire the elastic constant runs
elif [ $test_type == elastic ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        ElasticPrepFire.sh $cryst_sys
        cd ../..
    done

# display the fitting results and solve the elastic constants
elif [ $test_type == solve ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        ElasticDispSolveMod.sh $cryst_sys
        cd ../..
    done

# prep and fire a different strain
elif [ $test_type == alternative ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    if [ $cryst_sys == cubic ]; then
        alternative='2c11+2c12+c44'
    fi
    for n in $element_list
    do
        cd $n/$cryst_struct/elastic || exit 1
        Prepare.sh $alternative $cryst_sys
        Fire.sh $alternative
        cd ../../..
    done

# solve the elastic constants with the new set of strains
elif [ $test_type == solve_again ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    if [ $cryst_sys == cubic ]; then
        alternative='2c11+2c12+c44'
    fi
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        (cd elastic; Display.sh $alternative)
        ElasticSolve.sh $cryst_sys alternative
        cd ../..
    done

fi
