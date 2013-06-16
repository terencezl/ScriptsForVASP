#!/bin/bash
# Use: In the top working directory
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# Elastic.sh cubic / tetragonal /... prep-fire / disp-solve

if [ $1 == cubic ]; then
    dir_list="c11+2c12 c11-c12 c44"
elif [ $1 == cubic_A ]; then
    dir_list="A1 A2 A3 A4 A5 A6"
elif [ $1 == tetragonal ]; then echo
    dir_list="c11 c33 c44 5c11-4c12-2c13+c33 c11+c12-4c13+2c33 c11+c12-4c13+2c33+2c66"
elif [ $1 == orthorhombic ]; then echo
    dir_list="c11 c22 c33 c44 c55 c66 4c11-4c12-4c13+c22+2c23+c33 c11-4c12+2c13+4c22-4c23+c33 c11+2c12-4c13+c22-4c23+4c33"
elif [ $1 == hexagonal ]; then echo
elif [ $1 == trigonal ]; then echo
elif [ $1 == monoclinic ]; then echo
elif [ $1 == triclinic ]; then echo
fi

if [ $2 == prep-fire ]; then
    mkdir elastic 2> /dev/null
    cd elastic
    cp -r ../INPUT .

    scaling_factor=$(grep "Equilibrium scaling factor is" ../lctest/lctest_output.txt | head -1 | awk '{print $5}')
    sed -i "s/@R@/$scaling_factor/g" INPUT/POSCAR
    sed -i 's/#NSW /NSW /g' INPUT/INCAR
    sed -i 's/#IBRION /IBRION /g' INPUT/INCAR
    sed -i 's/#EDIFFG /EDIFFG /g' INPUT/INCAR
    #sed -i 's/ISMEAR /#ISMEAR /g' INPUT/INCAR

    for n in $dir_list
    do
        Prep-fire.sh $n $1
    done

elif [ $2 == disp-solve ]; then
    cd elastic || exit 1

    for n in $dir_list
    do
        Display.sh $n
    done

    echo $PWD | tee elastic_output.txt
    cd ..
    Elastic-solve.sh $1 original
    
fi
