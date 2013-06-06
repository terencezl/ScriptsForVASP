#!/bin/bash
# Use: In the top working directory,
# Display.sh TEST_TYPE (e.g. entest/c11-c12) P (polynomial fit)/M (Murnaghan fit)

if [[ "$1" == */ ]]; then test_type=${1%/}; else test_type=$1; fi
cd $test_type || exit 1
fname=$test_type"_output.txt"
echo $PWD > $fname
for n in $(ls -F)
do
    if [[ "$n" == */ ]]; then dir_list=$dir_list" "${n%/}; fi
done

data_line_count=0
if [[ $test_type == "entest" || $test_type == "kptest" ]]; then
    for n in $dir_list
    do
        if [[ "$n" != *-1 ]]; then dir_list_enkp=$dir_list_enkp" "$n; fi
    done
    dir_list=$(echo -e ${dir_list_enkp// /\\n} | sort -n)
    
    smallest=$(echo $dir_list | awk '{print $1}')
    largest=$(echo $dir_list | awk '{print $NF}')
    step=$(( $(echo $dir_list | awk '{print $2}') - $smallest ))
    lc1=$(sed -n '2p' $smallest/POSCAR)
    lc2=$(echo "$lc1+0.1" | bc)
    
    if [ $test_type == "entest" ]; then
        echo -e "ENCUT from $smallest to $largest step $step" >> $fname
        echo -e "LC1 = $lc1, LC2 = $lc2 (directories with \"-1\")" >> $fname
        echo -e "\nENCUT(eV)\tE_LC1(eV)\tE_LC2(eV)\tDE=E_LC2-E_LC1(eV)\tdDE(eV)" >> $fname
    else
        echo -e "nKP from $smallest to $largest step $step" >> $fname
        echo -e "LC1 = $lc1, LC2 = $lc2 (directories with \"-1\")" >> $fname
        echo -e "\nnKP\tE_LC1(eV)\tE_LC2(eV)\tDE=E_LC2-E_LC1(eV)\tdDE(eV)" >> $fname
    fi
    
    for n in $dir_list
    do
        E_LC1=$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')        # get the total energy when LC1
        E_LC2=$(grep sigma $n-1/OUTCAR | tail -1 | awk '{print $7;}')      # get the total energy when LC2, characterized by patterns like 160-1
        DE=$(echo "($E_LC2)-($E_LC1)" | bc)                                # echo the expr into bc (calculator) to get DE. bash doesn't support floats
        if [ $n == $smallest ]; then DE_pre=$DE; fi                        # in the first cycle let dDE = 0
        dDE=$(echo "($DE)-($DE_pre)" | bc)                                 # get dDE (the diff of energy difference between two adjacent rows)
        dDE=${dDE#-}                                                       # return the absolute value of dDE
        echo -e "$n\t$E_LC1\t$E_LC2\t$DE\t$dDE" >> $fname                  # output into a file
        DE_pre=$DE                                                         # set DE for this cycle as DE_pre for the next cycle to get the next dDE
        data_line_count=$(($data_line_count + 1))
    done
    
#    Another routine
#    x=$(echo $(sed -n 8,$((7 + $data_line_count))p $fname |awk '{print $1}'))
#    x=[$(echo ${x// /,})]
    _Display_fit.py $test_type 5 $((5 + data_line_count)) >> $fname

elif [[ $test_type == "lctest" || $test_type == "rttest" || $test_type == "mesh2d" ]]; then
    dir_list=$(echo -e ${dir_list// /\\n} | sort -n)
    smallest=$(echo $dir_list | awk '{print $1}')
    largest=$(echo $dir_list | awk '{print $NF}')
    step=$(echo "$(echo $dir_list | awk '{print $2}') - $smallest" | bc)
    if [[ $step == .* ]]; then step=0$step; fi

    if [ $test_type == "lctest" ]; then
        echo -e "Scaling factor from $smallest to $largest step $step" >> $fname
        echo -e "\nScalingFactor(Ang)\tVolume(Ang^3)\tE(eV)" >> $fname
    elif [ $test_type == "rttest" ]; then
        echo -e "Ratio from $smallest to $largest step $step" >> $fname
        echo -e "\nRatio\tE(eV)" >> $fname
#    needs more work
#    elif [ $test_type == "mesh2d" ]
#        echo "Scaling factor from $2 to $3 step $6, ratio from $4 to $5 step $6" >> $fname
#        echo -e "\nScalingFactor(Ang)\tRatio\tE(eV)" >> $fname
    fi
    
    for n in $dir_list
    do
        Vpcell=$(cat $n/OUTCAR |grep 'volume of cell' | tail -1 | awk '{print $5;}')
        echo -e "$n\t$Vpcell\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_rms=$(grep "FORCES:" $n/OUTCAR | tail -1 | awk '{print $6}')
        force_rms=${force_rms#-}
        if [ $(echo "$force_rms > 0.01" | bc) == 1 ]; then force_converge_flag=$n; fi
        entropy=$(grep "entropy T\*S" $n/OUTCAR | tail -1 | awk '{print $5}')
        entropy=${entropy#-}
        if [ $(echo "$entropy > 0.001" | bc) == 1 ]; then entropy_converge_flag=$n; fi
    done
    
    echo '' >> $fname
    _Display_fit.py $test_type 4 $((4 + data_line_count)) $2 >> $fname
    
    echo $PWD
    if [ $force_converge_flag ]; then echo "!Force doesn't converge during $force_converge_flag run!" | tee -a $fname; fi
    if [ $entropy_converge_flag ]; then echo "!Entropy doesn't converge during $entropy_converge_flag run!" | tee -a $fname; fi
    grep "R-squared is" $fname
    grep "Equilibrium scaling factor is" $fname
    grep "!Equilibrium point is out of the considered range!" $fname
    grep "B0 =" $fname
    grep "B0' =" $fname
    grep "Total energy is" $fname

elif [[ $test_type == *c[1-9][1-9]* ]]; then                              # meaning elastic const.
    echo "Delta from -0.03 to 0.03 with step 0.01" >> $fname
    echo -e "\nDelta(ratio)\tE(eV)" >> $fname
    dir_list="0.03n 0.02n 0.01n 0.00 0.01 0.02 0.03"
    
    for n in $dir_list
    do
        if [[ "$n" == *n ]]; then i=-${n%n}; else i=$n; fi
        echo -e "$i\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_rms=$(grep "FORCES:" $n/OUTCAR | tail -1 | awk '{print $6}')
        force_rms=${force_rms#-}
        if [ $(echo "$force_rms > 0.01" | bc) == 1 ]; then force_converge_flag=$i; fi
        entropy=$(grep "entropy T\*S" $n/OUTCAR | tail -1 | awk '{print $5}')
        entropy=${entropy#-}
        if [ $(echo "$entropy > 0.001" | bc) == 1 ]; then entropy_converge_flag=$i; fi
    done
    echo '' >> $fname
    _Display_fit.py $test_type 4 $((4 + data_line_count)) >> $fname
    echo $PWD
    if [ $force_converge_flag ]; then echo "!Force doesn't converge during $force_converge_flag run!" | tee -a $fname; fi
    if [ $entropy_converge_flag ]; then echo "!Entropy doesn't converge during $entropy_converge_flag run!" | tee -a $fname; fi
    grep "R-squared is" $fname
    
fi

echo -e "\nMaximal time per run: \c" >> $fname
largest=$(echo $dir_list | awk '{print $NF}')
cat < $(find $largest -mindepth 1 -name "*$largest*") | grep real | awk '{print $2;}' >> $fname
