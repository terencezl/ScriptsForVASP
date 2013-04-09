#!/bin/bash
# All batch files and the four input files should be prepared in the top working directory.
# In the top working directory:
# ./prepare.sh entest min max step nKP LC
# ./prepare.sh kptest min max step ENCUT LC
# ./prepare.sh lctest min max step ENCUT nKP
# Change the value of ENCUT of INCAR, nKP of KPOINTS and generalized length of POSCAR to @R@ before executing the create batch files.

if [[ $1 == "entest" || $1 == "kptest" ]]; then
    mkdir $1
    cd $1
    for ((n=$2; n<=$3; n=n+$4))
    do
        for i in $n $n-1
        do
            mkdir $i
            cd $i
            cp ../../INCAR .
            cp ../../POSCAR .
            cp ../../POTCAR .
            cp ../../KPOINTS .
            cp ../../qsub.parallel .
            if [ $1 == "entest" ]; then
                sed -i s/@R@/$n/g INCAR
    	        sed -i s/@R@/$5/g KPOINTS
            else
                sed -i s/@R@/$n/g KPOINTS
                sed -i s/@R@/$5/g INCAR
            fi
            sed -i s/@R@/$6/g POSCAR
            PWD=$(pwd)
            sed -i s%@R@%$PWD%g qsub.parallel
            cd ..
        done
    done
elif [ $1 == "lctest" ]; then
    mkdir lctest
    cd lctest
    List=""
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")
    do i=$(echo "scale=2;$n/1" | bc)			# change decimal format from 5.1 to 5.10
    List=$List" "$i
    done
    for n in $List
    do
        mkdir $n
        cd $n
        cp ../../INCAR .
        cp ../../POSCAR .
        cp ../../POTCAR .
        cp ../../KPOINTS .
        cp ../../qsub.parallel .
        sed -i s/@R@/$n/g POSCAR
	sed -i s/@R@/$5/g INCAR
	sed -i s/@R@/$6/g KPOINTS
        PWD=$(pwd)
        sed -i s%@R@%$PWD%g qsub.parallel
        cd ..
    done
else
    echo "Specify what you are going to test! entest/kptest/lctest"
fi
