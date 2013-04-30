#!/bin/bash
# Use: In the top working directory,
# ./Fire.sh TEST_TYPE min max step

if [[ $1 == "entest" || $1 == "kptest" ]]; then
    cd $1                                                    # get into the test dir
    for ((n=$2; n<=$3; n=n+$4))                              # submit each subfolder's command
    do
        (cd $n; qsub qsub.parallel)
        (cd $n-1; qsub qsub.parallel)
    done
elif [ $1 == "lctest" ]; then
    cd $1                                                    # get into the test dir
#    List=""                                                  # clear the initial variable
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")    # generate subfolders specified by float numbers
    do
        i=$(echo "scale=2;$n/1" | bc)                        # change decimal format from 5.1 to 5.10
        List=$List" "$i                                      # add List of float numbers up
    done
    for n in $List
    do                                                       # submit each subfolder's command
        (cd $n; qsub qsub.parallel)
    done
elif [ $1 == "rttest" ]; then
    cd $1                                                    # get into the test dir
#    List=""                                                  # clear the initial variable
    for n in $(awk "BEGIN{for(i=$2;i<=$3;i+=$4)print i}")    # generate subfolders specified by float numbers
    do
        i=$(echo "scale=2;$n/1" | bc)                        # change decimal format from 5.1 to 5.10
        List=$List" "$i                                      # add List of float numbers up
    done
    for n in $List
    do                                                       # submit each subfolder's command
        (cd $n; qsub qsub.parallel)
    done
else
    echo "Specify what you are going to test! entest/kptest/lctest"
fi
