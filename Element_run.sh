#!/usr/bin/env bash
# Use: In the TMN working directory where there is a list of folders by the name of elements
# Change element.dat first

element_list=$(echo $(cat $1))
echo "Tell me what you wish to do: prepare / lctest / display / equi-relax / copy-CONTCAR / elastic / solve / scrun / dosrun / plot-ldos / plot-tdos / bader-prerun / bader"
read -r test_type

# create INPUT folder in the selected cryst structure of each element
if [ $test_type == prepare ]; then
    for n in $element_list
    do
        mkdir -p "$n"/INPUT
        cp -r INPUT_ELEMENT/{INCAR,KPOINTS,POSCAR,qsub.parallel,qbader.serial} "$n"/INPUT
        cat ~/terencelz/local/potential-database/PAW-GGA/POTCAR_{"$n",N} > "$n"/INPUT/POTCAR
        cd "$n"/INPUT
        sed -i "/SYSTEM/c SYSTEM = $n" INCAR
        sed -i "1c $n" POSCAR
        cd ../..
    done

# basic lctest in the "top working directory"
elif [ $test_type == lctest ]; then
    echo "Enter the starting scaling factor, ending scaling factor and the step, separated by a space."
    read -r start end step
    for n in $element_list
    do
        cd "$n" || exit 1
        Prep-fire.sh lctest $start $end $step
        cd ..
    done

# fit the lctest curve to Murnaghan eqn of state
elif [ $test_type == display ]; then
    for n in $element_list
    do
        cd "$n" || exit 1
        Display.sh lctest
        echo
        cd ..
    done

# create an equilibrium run
elif [ $test_type == equi-relax ]; then
    for n in $element_list
    do
        cd "$n" || exit 1
        Fast_prep.sh equi-relax
        F equi-relax
        cd ..
    done

elif [ $test_type == copy-CONTCAR ];then
    for n in $element_list
    do
        cd "$n" || exit 1
        Display.sh equi-relax
        cd ..
    done

# prep and fire the elastic constant runs
elif [ $test_type == elastic ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd "$n" || exit 1
        Elastic.sh prep-fire $cryst_sys
        cd ..
    done

# display the fitting results and solve the elastic constants
elif [ $test_type == solve ]; then
    echo "Tell me the crystallographic system the compound is in: cubic / tetragonal / orthorhombic/..."
    read -r cryst_sys
    for n in $element_list
    do
        cd "$n" || exit 1
        Elastic.sh disp-solve $cryst_sys
        echo
        cd ..
    done

# do the self-consistent run, find k-point dosrun, bandstructure run, and plot them."
elif [ $test_type == scrun ] || [ $test_type == dosrun ] || [ $test_type == bsrun ] || [ $test_type == plot-ldos ] || [ $test_type == plot-tdos ]; then
    if [ $test_type == plot-ldos ] || [ $test_type == plot-tdos ]; then
        echo "Tell me the size of the frame. e.g. [-20,20,-5,5] for ldos, [-20,20,0,10] for tdos."
        read -r frame_size
    fi
    for n in $element_list
    do
        cd "$n" || exit 1
        Electronic.sh $test_type "$frame_size"
        cd ..
    done

elif [ $test_type == bader-prerun ]; then
    for n in $element_list
    do
        cd "$n"
        Bader.sh prerun
        cd ..
    done
    
elif [ $test_type == bader ]; then
    for n in $element_list
    do
        cd "$n"
        Bader.sh bader
        cd ..
    done

fi
