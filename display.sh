#!/bin/bash
# Use: In the top working directory,
# Display.sh TEST_TYPE (e.g. entest/c11-c12)

if [[ "$1" == */ ]]; then test_type=${1%/}; else test_type=$1; fi
cd $test_type || exit 1
fname=$test_type"_output.txt"
for n in $(ls -F)
do
    if [[ "$n" == */ ]]; then dir_list=$dir_list" "${n%/}; fi
done

data_line_count=0
if [[ $test_type == "entest" || $test_type == "kptest" ]]; then
    if [ $test_type == "entest" ]; then echo -e "\nENCUT(meV)\tE_LC1\tE_LC2\tDE=E_LC2-E_LC1\tdDE" >> $fname
    else echo -e "\nnKP\tE_LC1\tE_LC2\tDE=E_LC2-E_LC1\tdDE" >> $fname
    fi
    for n in $dir_list
    do
        if [[ "$n" != *-1 ]]; then dir_list_enkp=$dir_list_enkp" "$n; fi
    done
    dir_list=$(echo -e ${dir_list_enkp// /\\n} | sort -n)
    smallest=$(echo $dir_list | awk '{print $1}')
    for n in $dir_list
    do
        E_LC1=$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')        # get the total energy when LC1
        E_LC2=$(grep sigma $n-1/OUTCAR | tail -1 | awk '{print $7;}')      # get the total energy when LC2, characterized by patterns like 160-1
        DE=$(echo "($E_LC2)-($E_LC1)" | bc)                         # echo the expr into bc (calculator) to get DE. bash doesn't support floats
        if [ $n == $smallest ]; then DE_pre=$DE; fi                        # in the first cycle let dDE = 0
        dDE=$(echo "($DE)-($DE_pre)" | bc)                          # get dDE (the diff of energy difference between two adjacent rows)
        dDE=${dDE#-}                                                # return the absolute value of dDE
        echo -e "$n\t$E_LC1\t$E_LC2\t$DE\t$dDE" >> $fname          # output into a file
        DE_pre=$DE                                                  # set DE for this cycle as DE_pre for the next cycle to get the next dDE
        data_line_count=$(($data_line_count + 1))
    done
#    Another routine
#    x=$(echo $(sed -n 8,$((7 + $data_line_count))p $fname |awk '{print $1}'))
#    x=[$(echo ${x// /,})]
    _Display_fit.py $test_type 5 $((5+data_line_count)) >> $fname

elif [[ $test_type == "lctest" || $test_type == "rttest" || $test_type == "mesh2d" ]]; then
    if [ $test_type == "lctest" ]; then echo -e "\nScalingFactor (Ang)\tVolume (Ang^3)\tE" >> $fname
    elif [ $test_type == "rttest" ]; then echo -e "\nRatio\tE" >> $fname
    else echo -e "\nScalingFactor (Ang)\tRatio\tE" >> $fname
    fi
    for n in $dir_list
    do
        Vpcell=$(cat $n/OUTCAR |grep 'volume of cell' |tail -1| awk '{print $5;}')
        echo -e "$n\t$Vpcell\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
    done
    echo '' >> $fname
    _Display_fit.py $test_type 4 $((4+data_line_count)) >> $fname

elif [[ $test_type == *c[1-9][1-9]* ]]; then                              # meaning elastic const.
    echo -e "\nDelta (ratio)\tE" >> $fname
    dir_list="0.04n 0.03n 0.02n 0.01n 0.00 0.01 0.02 0.03 0.04"
    for n in $dir_list
    do
        if [[ "$n" == *n ]]; then i=-${n%n}; else i=$n; fi
        echo -e "$i\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
    done
    echo '' >> $fname
    _Display_fit.py $test_type 5 $((5+data_line_count)) >> $fname
    
fi

echo -e "Maximal time per run: \c" >> $fname
largest=$(echo $dir_list | awk '{print $NF}')
cat < $(find $largest -mindepth 1 -name "*$largest*") | grep real | awk '{print $2;}' >> $fname
