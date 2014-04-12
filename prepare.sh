#!/bin/bash
# The four input files and qsub.parallel should be prepared in the top working directory.
# In the top working directory:
# Prepare.sh entest min max step LC1
# Prepare.sh kptest min max step LC1
# Prepare.sh lctest min max step
# Prepare.sh rttest min max step
# Prepare.sh mesh2d LC_min LC_max RT_min RT_max step
# Prepare.sh c11-c12 cubic
# Note: Change the value of ENCUT of INCAR (when doing entest), nKP of KPOINTS (when kptest) and scaling factor of POSCAR (when both) to @R@ before executing the batch file.

function qsub_replacer {
    #qname=${PWD//\//_}
    #qname=${qname##*utl0268_}
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
    sed -i s/@N@/$qname/g qsub.parallel
    sed -i s%@R@%$PWD%g qsub.parallel
}

test_type=$1
mkdir "$test_type" 2> /dev/null
cd "$test_type" || exit 1
test_type=${test_type%%_*}
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
            qsub_replacer
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
        qsub_replacer
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
        qsub_replacer
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
        qsub_replacer
        cd ..
    done

#elif [[ "$test_type" == "equi-relax" ]]; then
#    cd ..
#    scaling_factor=$(grep "Equilibrium scaling factor is" lctest/lctest_output.txt | head -1 | awk '{print $5}')
#    Prep-fast.sh $test_type
#    sed -i "2c $scaling_factor" equi-relax/POSCAR

#elif [ "$test_type" == "mesh2d" ]; then
#    for x in $(awk "BEGIN{for(i=$2;i<=$3;i+=$6)print i}")                   # generate subfolders specified by float numbers
#    do i=$(echo "scale=2;$x/1" | bc)                                        # change decimal format from 5.1 to 5.10
#    dir_listx=$dir_listx" "$i                                                         # add dir_list of float numbers up
#    done
#    for y in $(awk "BEGIN{for(i=$4;i<=$5;i+=$6)print i}")                   # generate subfolders specified by float numbers
#    do i=$(echo "scale=2;$y/1" | bc)                                        # change decimal format from 5.1 to 5.10
#    dir_listy=$dir_listy" "$i                                                         # add dir_list of float numbers up
#    done
#    for x in $dir_listx
#    do
#        for y in $dir_listy
#        do
#            mkdir x$x-y$y
#            cd x$x-y$y
#            cp ../../INPUT/INCAR .
#            cp ../../INPUT/POSCAR .
#            cp ../../INPUT/POTCAR .
#            cp ../../INPUT/KPOINTS .
#            cp ../../INPUT/qsub.parallel .
#            sed -i s/@Rx@/$x/g POSCAR                                            # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS/POSCAR are reserved
#            sed -i s/@Ry@/$y/g POSCAR                                            # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS/POSCAR are reserved
#            qsub_replacer
#            cd ..
#        done
#    done
    
elif [[ $test_type == *c[1-9][1-9]* || $test_type == A* ]]; then
    if [ $test_type == c44 ]; then
#        dir_list="0.040n 0.025n 0.025 0.040"
#        dir_list="0.050n 0.030n 0.030 0.050"
        dir_list="0.030 0.050n"
    elif [ $test_type == c11-c12 ]; then
        dir_list="0.025 0.040n"
#        dir_list="0.015"
#        dir_list="0.040n 0.025n 0.025 0.040"
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
        qsub_replacer
        if [[ "$n" == *n ]]; then n=-${n%n}; fi
        _Prepare-strain.py $test_type $2 $n
        cd ..
    done

else
    echo "Specify what you are going to test!"
fi
