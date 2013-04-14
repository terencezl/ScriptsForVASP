#!/bin/bash
# All batch files and the four input files should be prepared in the top working directory.
# In the top working directory:
# ./prepare.sh entest min max step nKP LC1
# ./prepare.sh kptest min max step ENCUT LC1
# ./prepare.sh lctest min max step
# ./prepare.sh rttest min max step
# Change the value of ENCUT of INCAR, nKP of KPOINTS and generalized length of POSCAR to @R@ before executing the create batch files.

if [[ $1 == "entest" || $1 == "kptest" ]]; then
#   ls $1 || mkdir $1 && cd $1
    mkdir $1 2> /dev/null                                                   # if wishing to add more trials for a test but test dir already exists, ignore this error
    cd $1                                                                   # get into the test dir
    fname=$1"_output.txt"
    lc=$(echo $6+0.1 | bc)                                                  # get LC2=LC1+0.1. bc is calculator; bash doesn't support floats
    echo -e $1\\n > $fname                                                  # start to write some head info of each trial run
    if [ $1 == "entest" ]; then
        echo -e ENCUT from $2 to $3 step $4\\nnKP = $5 > $fname
    else
        echo -e nKP from $2 to $3 step $4\\nENCUT = $5 >> $fname
    fi
    echo -e "LC1 = $6, LC2 = $lc (directories with \"-1\")" >> $fname
    for ((n=$2; n<=$3; n=n+$4))                                             # create each subfolder
    do
        for i in $n $n-1                                                    # subfolders for two LCs
        do
            mkdir $i
            cd $i
            cp ../../INCAR .
            cp ../../POSCAR .
            cp ../../POTCAR .
            cp ../../KPOINTS .
            cp ../../qsub.parallel .
            if [ $1 == "entest" ]; then
                sed -i s/@R@/$n/g INCAR                                     # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS are reserved
                sed -i s/@R@/$5/g KPOINTS
            else
                sed -i s/@R@/$n/g KPOINTS
                sed -i s/@R@/$5/g INCAR
            fi
            if [ $i == $n ]; then                                           # replace the two LCs in their own POSCAR. Arguments about the generalized length is reserved
                sed -i s/@R@/$6/g POSCAR
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
elif [ $1 == "lctest" ]; then
    mkdir $1 2> /dev/null
    cd $1                                                                   # get into the test dir
    fname=$1"_output.txt"
    echo -e $1\\n > $fname                                                  # start to write some head info of each trial run
    echo LC from $2 to $3 step $4 >> $fname
#    echo -e ENCUT = $5\\nnKP = $6 >> $fname
    List=""                                                                 # clear the initial variable
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")                   # generate subfolders specified by float numbers
    do i=$(echo "scale=2;$n/1" | bc)                                        # change decimal format from 5.1 to 5.10
    List=$List" "$i                                                         # add List of float numbers up
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
        sed -i s/@R@/$n/g POSCAR                                            # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS/POSCAR are reserved
 #       sed -i s/@R@/$5/g INCAR
 #       sed -i s/@R@/$6/g KPOINTS
        qname=${PWD//\//_}                                                  # edit the name of the command so that it looks like terencelz_GaN_MN_lctest_140
        qname=${qname##*utl0268_}
        sed -i s/@N@/$qname/g qsub.parallel                                 # replace the name in the command file. Arg reserved
        sed -i s%@R@%$PWD%g qsub.parallel                                   # replace the trial subfolder in the command file. Arg reserved
        cd ..
    done
elif [ $1 == "rttest"]; then
    mkdir $1 2> /dev/null
    cd $1                                                                   # get into the test dir
    fname=$1"_output.txt"
    echo -e $1\\n > $fname                                                  # start to write some head info of each trial run
    echo "Ratio c/a from $2 to $3 step $4" >> $fname
#    echo -e ENCUT = $5\\nnKP = $6 >> $fname
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")                   # generate subfolders specified by float numbers
    do
        i=$(echo "scale=2;$n/1" | bc)                                       # change decimal format from 5.1 to 5.10
        List=$List" "$i                                                     # add List of float numbers up
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
        sed -i s/@RT@/$n/g POSCAR                                            # use sed -i s/xx/yy/g FILE to do replacement. Arguments in INCAR/KPOINTS/POSCAR are reserved
#        sed -i s/@R@/$5/g INCAR
#        sed -i s/@R@/$6/g KPOINTS
        qname=${PWD//\//_}                                                  # edit the name of the command so that it looks like terencelz_GaN_MN_lctest_140
        qname=${qname##*utl0268_}
        sed -i s/@N@/$qname/g qsub.parallel                                 # replace the name in the command file. Arg reserved
        sed -i s%@R@%$PWD%g qsub.parallel                                   # replace the trial subfolder in the command file. Arg reserved
        cd ..
    done
else
    echo "Specify what you are going to test! entest/kptest/lctest"
fi
