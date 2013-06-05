#!/bin/bash
# Use: In the TMN working directory
# Change element.dat first

element_list=$(cat INPUT_ELEMENT/element.dat)
echo "Tell me what you wish to do: prepare/lctest/display/elastic/solve/alternative/solve_again"
read -r test_type
echo "What lattice type are you thinking of? zinc-blende/rocksalt/CsCl/..."
read -r cryst_struct

if [ $test_type == prepare ]; then
    for n in $element_list
    do
        mkdir -p $n/$cryst_struct/INPUT
        cp -r INPUT_ELEMENT/{INCAR,KPOINTS,POSCAR_$cryst_struct,qsub.parallel} $n/$cryst_struct/INPUT
        cat INPUT_ELEMENT/POTCAR_LDA/{POTCAR_$n,POTCAR_N} > $n/$cryst_struct/INPUT/POTCAR
        cd $n/$cryst_struct/INPUT
        mv POSCAR_$cryst_struct POSCAR
        sed -i "s/@N@/$n N ($cryst_struct)/g" INCAR
        sed -i "s/@N@/$n N ($cryst_struct)/g" POSCAR
        cd ../../..
    done

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

elif [ $test_type == display ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        Display.sh lctest
        cd ../..
    done

elif [ $test_type == elastic ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        ElasticPrepFire.sh $cryst_sys
        cd ../..
    done

elif [ $test_type == solve ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct || exit 1
        ElasticDispSolve.sh $cryst_sys
        cd ../..
    done

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
