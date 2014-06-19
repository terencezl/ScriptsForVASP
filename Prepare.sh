#!/bin/bash
# The four input files and qsub.parallel should be prepared in the top working directory.
# Prepare.sh entest min max step LC1
# Prepare.sh kptest min max step LC1
# Prepare.sh lctest min max step
# Prepare.sh rttest min max step
# Prepare.sh c11-c12 cubic

function qsub_replacer {
    qname_1=${PWD##*/}
    PWD_2=${PWD%/*}
    qname_2=${PWD_2##*/}
    PWD_3=${PWD_2%/*}
    qname_3=${PWD_3##*/}
    PWD_4=${PWD_3%/*}
    qname_4=${PWD_4##*/}
    qname="T$qname_4$qname_3$qname_2$qname_1"
    if [[ $(echo $qname | wc -c) > 17 ]]; then
        qname="T$qname_3$qname_2$qname_1"
    fi
    sed -i "/#PBS -N/c #PBS -N $qname" $1
    sed -i "/^cd/c cd $PWD" $1
}

directory_name=$1
mkdir "$directory_name" 2> /dev/null
cd "$directory_name" || exit 1
test_type=${directory_name%%_*}
fname="$test_type""_output.txt"

if [[ "$test_type" == "entest" || "$test_type" == "kptest" ]]; then
    lc1=$5
    lc2=$(echo "$lc1+0.1" | bc)                                                  # get LC2=LC1+0.1. bc is calculator; bash doesn't support floats
    for ((n=$2; n<=$3; n=n+$4))                                             # create each subfolder
    do
        for i in $n $n-1                                                    # subfolders for two LCs
        do
            mkdir $i
            cd $i
            cp ../../INPUT/INCAR .
            cp ../../INPUT/POSCAR .
            cp ../../INPUT/POTCAR .
            cp ../../INPUT/KPOINTS .
            cp ../../INPUT/qsub.parallel .
            if [ "$test_type" == "entest" ]; then
                sed -i "s/.*ENCUT.*/ENCUT = $n/g" INCAR
            else
                sed -i "4c $n $n $n" KPOINTS
            fi
            if [ $i == $n ]; then                                           # replace the two LCs in their own POSCAR. Arguments about the generalized length is reserved
                sed -i "2c $lc1" POSCAR
            else
                sed -i "2c $lc2" POSCAR
            fi
            qsub_replacer qsub.parallel
            cd ..
        done
    done

elif [[ "$test_type" == "lctest" ]]; then
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")                   # generate subfolders specified by float numbers
    do
        i=$(echo "scale=3;$n/1" | bc)                                        # change decimal format from 5.1 to 5.10
        dir_list=$dir_list" "$i                                                         # add dir_list of float numbers up
    done
    for n in $dir_list
    do
        mkdir $n
        cd $n
        cp ../../INPUT/INCAR .
        cp ../../INPUT/POSCAR .
        cp ../../INPUT/POTCAR .
        cp ../../INPUT/KPOINTS .
#        cp ../../INPUT/WAVECAR .
        cp ../../INPUT/qsub.parallel .
        sed -i "2c $n" POSCAR
        qsub_replacer qsub.parallel
        cd ..
    done

elif [[ "$test_type" == "rttest" ]]; then
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")                   # generate subfolders specified by float numbers
    do
        i=$(echo "print('{0:.3f}'.format($n))" | python)                              # change decimal format from 5.1 to 5.10
        dir_list=$dir_list" "$i                                                         # add dir_list of float numbers up
    done
    for n in $dir_list
    do
        mkdir $n
        cd $n
        cp ../../INPUT/INCAR .
        cp ../../INPUT/POSCAR .
        cp ../../INPUT/POTCAR .
        cp ../../INPUT/KPOINTS .
        cp ../../INPUT/WAVECAR .
        cp ../../INPUT/qsub.parallel .
        sed -i "s/@R@/$n/g" POSCAR
        qsub_replacer qsub.parallel
        cd ..
    done

elif [[ "$test_type" == "agltest" ]]; then
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")
    do
        n=$(echo "print('{0:f}'.format($n))" | python)
        if [[ "$n" == -* ]]; then n=${n#-}n; fi
        dir_list=$dir_list" "$n
    done
    for n in $dir_list
    do
        mkdir $n
        cd $n
        cp ../../INPUT/INCAR .
        cp ../../INPUT/POSCAR .
        cp ../../INPUT/POTCAR .
        cp ../../INPUT/KPOINTS .
        cp ../../INPUT/qsub.parallel .
        if [[ "$n" == *n ]]; then n=-${n%n}; fi
        _ions-position-rotator.py $n
        qsub_replacer qsub.parallel
        cd ..
    done


elif [[ $test_type == *c[1-9][1-9]* || $test_type == A* ]]; then
    if [ $test_type == c44 ]; then
        dir_list="0.020 0.035 0.050n"
    elif [ $test_type == c11-c12 ]; then
        dir_list="0.020 0.030 0.040n"
    fi
    for n in $dir_list
    do
        mkdir $n
        cd $n
        cp ../../INPUT/INCAR .
        cp ../../INPUT/POSCAR .
        cp ../../INPUT/POTCAR .
        cp ../../INPUT/KPOINTS .
        cp ../../INPUT/qsub.parallel .
        qsub_replacer qsub.parallel
        if [[ "$n" == *n ]]; then n=-${n%n}; fi
        _Prepare-strain.py $test_type $2 $n
        cd ..
    done

else
    echo "Specify what you are going to test!"
fi
