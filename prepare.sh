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

mkdir "$1" 2> /dev/null
cd "$1"
fname="$1""_output.txt"
echo -e "$1" >> $fname                                                  # start to write some head info of each trial run

if [[ "$1" == "entest" || "$1" == "kptest" ]]; then
    lc=$(echo $6+0.1 | bc)                                                  # get LC2=LC1+0.1. bc is calculator; bash doesn't support floats
    if [ "$1" == "entest" ]; then
        echo -e "ENCUT from $2 to $3 step $4" >> $fname
    else
        echo -e "nKP from $2 to $3 step $4" >> $fname
    fi
    echo -e "LC1 = $5, LC2 = $lc (directories with \"-1\")" >> $fname
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
            if [ "$1" == "entest" ]; then
                sed -i s/@R@/$n/g INCAR                                     # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS are reserved
            else
                sed -i s/@R@/$n/g KPOINTS
            fi
            if [ $i == $n ]; then                                           # replace the two LCs in their own POSCAR. Arguments about the generalized length is reserved
                sed -i s/@R@/$5/g POSCAR
            else
                sed -i s/@R@/$lc/g POSCAR
            fi
            qname=${PWD//\//_}                                              # edit the name of the command so that it looks like terencelz_GaN_MN_entest_140
            qname=${qname##*utl0268_}
            sed -i s/@N@/$qname/g qsub.parallel                             # replace the name in the command file. Arg reserved
            sed -i s%@R@%$PWD%g qsub.parallel                               # replace the trial subfolder in the command file. Arg reserved
            cd ..
        done
    done

elif [[ "$1" == "lctest" || "$1" == "rttest" ]]; then
    if [ "$1" == "lctest" ]; then
        echo -e "Scaling factor from $2 to $3 step $4" >> $fname
    else
        echo -e "Ratio from $2 to $3 step $4" >> $fname
    fi
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")                   # generate subfolders specified by float numbers
    do
        i=$(echo "scale=2;$n/1" | bc)                                        # change decimal format from 5.1 to 5.10
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
        cp ../../INPUT/qsub.parallel .
        sed -i s/@R@/$n/g POSCAR                                            # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS/POSCAR are reserved
#        sed -i s/@R@/$5/g INCAR
#        sed -i s/@R@/$6/g KPOINTS
        qname=${PWD//\//_}                                                  # edit the name of the command so that it looks like terencelz_GaN_MN_lctest_140
        qname=${qname##*utl0268_}
        sed -i s/@N@/$qname/g qsub.parallel                                 # replace the name in the command file. Arg reserved
        sed -i s%@R@%$PWD%g qsub.parallel                                   # replace the trial subfolder in the command file. Arg reserved
        cd ..
    done

elif [ "$1" == "mesh2d" ]; then
    echo "Scaling factor from $2 to $3 step $6, ratio from $4 to $5 step $6" >> $fname
    for x in $(awk "BEGIN{for(i=$2;i<=$3;i+=$6)print i}")                   # generate subfolders specified by float numbers
    do i=$(echo "scale=2;$x/1" | bc)                                        # change decimal format from 5.1 to 5.10
    dir_listx=$dir_listx" "$i                                                         # add dir_list of float numbers up
    done
    for y in $(awk "BEGIN{for(i=$4;i<=$5;i+=$6)print i}")                   # generate subfolders specified by float numbers
    do i=$(echo "scale=2;$y/1" | bc)                                        # change decimal format from 5.1 to 5.10
    dir_listy=$dir_listy" "$i                                                         # add dir_list of float numbers up
    done
    for x in $dir_listx
    do
        for y in $dir_listy
        do
            mkdir x$x-y$y
            cd x$x-y$y
            cp ../../INPUT/INCAR .
            cp ../../INPUT/POSCAR .
            cp ../../INPUT/POTCAR .
            cp ../../INPUT/KPOINTS .
            cp ../../INPUT/qsub.parallel .
            sed -i s/@Rx@/$x/g POSCAR                                            # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS/POSCAR are reserved
            sed -i s/@Ry@/$y/g POSCAR                                            # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS/POSCAR are reserved
            qname=${PWD//\//_}                                                  # edit the name of the command so that it looks like terencelz_GaN_MN_lctest_140
            qname=${qname##*utl0268_}
            sed -i s/@N@/$qname/g qsub.parallel                                 # replace the name in the command file. Arg reserved
            sed -i s%@R@%$PWD%g qsub.parallel                                   # replace the trial subfolder in the command file. Arg reserved
            cd ..
        done
    done
    
elif [[ $1 == *c[1-9][1-9]* ]]; then
    echo -e "Crystallographic system: $2" >> $fname
    echo "Delta from -0.04 to 0.04 with step 0.01" >> $fname
    dir_list="0.04n 0.03n 0.02n 0.01n 0.00 0.01 0.02 0.03 0.04"
    for n in $dir_list
    do
        mkdir $n
        cd $n
        cp ../../INPUT/INCAR .
        cp ../../INPUT/POSCAR .
        cp ../../INPUT/POTCAR .
        cp ../../INPUT/KPOINTS .
        cp ../../INPUT/qsub.parallel .
        qname=${PWD//\//_}
        qname=${qname##*utl0268_}
        sed -i s/@N@/$qname/g qsub.parallel
        sed -i s%@R@%$PWD%g qsub.parallel
        if [[ "$n" == *n ]]; then n=-${n%n}; fi
        _Prepare_strain.py $1 $2 $n
        cd ..
    done

else
    echo "Specify what you are going to test!"
fi
