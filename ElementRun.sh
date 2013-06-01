#!/bin/bash
# Use: In the TMN working directory
# Change element.dat first

element_list=$(cat INPUT/element.dat)
echo "Tell me what you wish to do: prepare/lctest/display/elastic/solve"
read -r test_type
echo "What lattice type are you thinking of? zinc-blende/rocksalt/..."
read -r cryst_struct

if [ $test_type == prepare ]; then
    for n in $element_list
    do
        mkdir -p $n/$cryst_struct/INPUT
        cp -r INPUT/INCAR INPUT/KPOINTS INPUT/POSCAR_$cryst_struct INPUT/qsub.parallel $n/$cryst_struct/INPUT
        cat INPUT/POTCAR_LDA/POTCAR_$n INPUT/POTCAR_LDA/POTCAR_N > $n/$cryst_struct/INPUT/POTCAR
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
        cd $n/$cryst_struct
        Prepare.sh lctest $start $end $step
        Fire.sh lctest $start $end $step
        cd ../..
    done

elif [ $test_type == display ]; then
    for n in $element_list
    do
        cd $n/$cryst_struct
        Display.sh lctest
        cd ../..
    done

elif [ $test_type == elastic ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct
        ElasticPrepFire.sh $cryst_sys
        cd ../..
    done

elif [ $test_type == solve ]; then
    echo "Tell me the crystallographic system the compound is in: cubic/tetragonal/orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd $n/$cryst_struct
        ElasticDispSolve.sh $cryst_sys
        cd ../..
    done

fi
