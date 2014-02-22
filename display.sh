#!/bin/bash
# Use: In the top working directory,
# Display.sh TEST_TYPE (e.g. entest/c11-c12)

function force_entropy_detector {
    force_max=$(grep "FORCES:" $1/OUTCAR | tail -1 | awk '{print $5}')
    force_max=${force_max#-}
    if [ $(echo "$force_max > 0.01" | bc) == 1 ]; then force_converge_flag="$force_converge_flag\n$1 $force_max"; fi
    entropy=$(grep "entropy T\*S" $1/OUTCAR | tail -1 | awk '{print $5}')
    entropy=${entropy#-}
    atom_sum_expr=$(echo $(sed -n '6p' $1/POSCAR))
    atom_sum=$((${atom_sum_expr// /+}))
    if [ $(echo "$entropy > 0.001 * $atom_sum" | bc) == 1 ]; then entropy_converge_flag="$entropy_converge_flag\n$1 $entropy"; fi
}

if [[ "$1" == */ ]]; then test_type=${1%/}; else test_type=$1; fi
cd $test_type || exit 1
test_type=${test_type%%_*}
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
        echo -e "\nENCUT(eV)\tE_LC1(eV)\tdE_LC1 (eV)\tE_LC2(eV)\tdE_LC2 (eV)\tDE=E_LC2-E_LC1(eV)\t\tdDE(eV)" >> $fname
    else
        echo -e "nKP from $smallest to $largest step $step" >> $fname
        echo -e "LC1 = $lc1, LC2 = $lc2 (directories with \"-1\")" >> $fname
        echo -e "\nnKP\tE_LC1(eV)\tdE_LC1 (eV)\tE_LC2(eV)\tdE_LC2 (eV)\tDE=E_LC2-E_LC1(eV)\t\tdDE(eV)" >> $fname
    fi
    
    for n in $dir_list
    do
        E_LC1=$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')        # get the total energy when LC1
        E_LC2=$(grep sigma $n-1/OUTCAR | tail -1 | awk '{print $7;}')      # get the total energy when LC2, characterized by patterns like 160-1
        DE=$(echo "($E_LC2)-($E_LC1)" | bc)                                # echo the expr into bc (calculator) to get DE. bash doesn't support floats
        if [ $n == $smallest ]; then
            DE_pre=$DE
            E_LC1_pre=$E_LC1
            E_LC2_pre=$E_LC2
        fi                                                                 # in the first cycle let dDE = 0
        dDE=$(echo "($DE)-($DE_pre)" | bc)                                 # get dDE (the diff of energy difference between two adjacent rows)
        dE_LC1=$(echo "($E_LC1)-($E_LC1_pre)" | bc)
        dE_LC2=$(echo "($E_LC2)-($E_LC2_pre)" | bc)
        dDE=${dDE#-}                                                       # return the absolute value of dDE
        dE_LC1=${dE_LC1#-}                                                       # return the absolute value of dDE
        dE_LC2=${dE_LC2#-}                                                       # return the absolute value of dDE
        echo -e "$n\t$E_LC1\t$dE_LC1\t$E_LC2\t$dE_LC2\t$DE\t\t$dDE" >> $fname                  # output into a file
        DE_pre=$DE                                                         # set DE for this cycle as DE_pre for the next cycle to get the next dDE
        E_LC1_pre=$E_LC1
        E_LC2_pre=$E_LC2
        data_line_count=$(($data_line_count + 1))
    done
    
    _Display-fit.py $test_type 5 $((5 + data_line_count)) >> $fname

elif [[ $test_type == "lctest" ]]; then
    dir_list=$(echo -e ${dir_list// /\\n} | sort -n)
    smallest=$(echo $dir_list | awk '{print $1}')
    largest=$(echo $dir_list | awk '{print $NF}')
    step=$(echo "print($(echo $dir_list | awk '{print $2}') - ($smallest))" | python)
    echo -e "Scaling factor from $smallest to $largest step $step" >> $fname
    echo -e "\nScalingFactor(Ang)\tVolume(Ang^3)\tE(eV)" >> $fname
    
    for n in $dir_list
    do
        Vpcell=$(cat $n/OUTCAR | grep 'volume of cell' | tail -1 | awk '{print $5;}')
        echo -e "$n\t$Vpcell\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_detector $n
    done
    
    echo $PWD
    echo '' >> $fname
    if [ "$force_converge_flag" ]; then echo -e "!Force doesn't converge during$force_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    if [ "$entropy_converge_flag" ]; then echo -e "!Entropy doesn't converge during$entropy_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    _Display-fit.py $test_type 4 $((4 + data_line_count)) >> $fname
    grep "!Equilibrium point is out of the considered range!" $fname
    grep "R-squared is" $fname
    grep "Equilibrium scaling factor is" $fname
    grep "B0 =" $fname
    grep "B0' =" $fname
    grep "Total energy is" $fname

    cd ..
    scaling_factor=$(grep "Equilibrium scaling factor is" lctest/lctest_output.txt | head -1 | awk '{print $5}')
    closest_lc=$(echo $dir_list | awk '{print $1}')
    for n in $dir_list; do
        diffn=$(echo "$n - $scaling_factor" | bc)
        diffn=${diffn#-}
        difflc=$(echo "$closest_lc - $scaling_factor" | bc)
        difflc=${difflc#-}
        if [[ $(echo "$diffn < $difflc" | bc) == 1 ]]; then closest_lc=$n; fi
    done
    cp lctest/$closest_lc/CONTCAR INPUT/POSCAR
    sed -i "2c $scaling_factor" INPUT/POSCAR
    cd $test_type

elif [[ $test_type == "rttest" ]]; then
    dir_list=$(echo -e ${dir_list// /\\n} | sort -n)
    smallest=$(echo $dir_list | awk '{print $1}')
    largest=$(echo $dir_list | awk '{print $NF}')
    step=$(echo "print($(echo $dir_list | awk '{print $2}') - ($smallest))" | python)

    echo -e "Ratio from $smallest to $largest step $step" >> $fname
    echo -e "\nRatio\tVolume\tE(eV)" >> $fname

    for n in $dir_list
    do
        Vpcell=$(cat $n/OUTCAR | grep 'volume of cell' | tail -1 | awk '{print $5;}')
        echo -e "$n\t$Vpcell\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_detector $n
    done
    
    echo $PWD
    echo '' >> $fname
    if [ "$force_converge_flag" ]; then echo -e "!Force doesn't converge during$force_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    if [ "$entropy_converge_flag" ]; then echo -e "!Entropy doesn't converge during$entropy_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    _Display-fit.py $test_type 4 $((4 + data_line_count)) >> $fname
    grep "!Equilibrium point is out of the considered range!" $fname
    grep "R-squared is" $fname
    grep "Equilibrium ratio is" $fname
    grep "B0 =" $fname
    grep "B0' =" $fname
    grep "Total energy is" $fname

elif [[ $test_type == "agltest" ]]; then
    dir_list=$(echo -e ${dir_list// /\\n} | sort -n)
    smallest=$(echo $dir_list | awk '{print $1}')
    largest=$(echo $dir_list | awk '{print $NF}')
    step=$(echo "print($(echo $dir_list | awk '{print $2}') - ($smallest))" | python)

    echo -e "Angle from $smallest to $largest step $step" >> $fname
    echo -e "\nAngle\tE(eV)" >> $fname

    for n in $dir_list
    do
        echo -e "$n\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_detector $n
    done
    
    echo $PWD
    echo '' >> $fname
    if [ "$force_converge_flag" ]; then echo -e "!Force doesn't converge during$force_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    if [ "$entropy_converge_flag" ]; then echo -e "!Entropy doesn't converge during$entropy_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    _Display-fit.py $test_type 4 $((4 + data_line_count)) >> $fname
    grep "R-squared is" $fname
    grep "Minimal total energy is" $fname

elif [[ $test_type == "equi-relax" ]]; then
    cp CONTCAR ../INPUT/POSCAR
    exit 0

elif [[ $test_type == *c[1-9][1-9]* || $test_type == A* ]]; then                              # meaning elastic const.
    if [[ $test_type == c44 ]]; then
        if [ -d 0.030n ]; then rm -r 0.030n; fi
        if [ -d 0.050 ]; then rm -r 0.050; fi
        cp -r 0.030 0.030n
        cp -r 0.050n 0.050
    elif [[ $test_type == c11-c12 ]]; then
        if [ -d 0.025n ]; then rm -r 0.025n; fi
        if [ -d 0.040 ]; then rm -r 0.040; fi
        cp -r 0.025 0.025n
#        cp -r 0.015 0.015n
        cp -r 0.040n 0.040
    fi
    if [ -d 0.000 ]; then rm -r 0.000; fi
    cp -r ../../equi-relax 0.000
    unset dir_list
    for n in $(ls -F)
    do
        if [[ "$n" == */ ]]; then dir_list=$dir_list" "${n%/}; fi
    done

    for n in $dir_list
    do
        if [[ "$n" == *n ]]; then i=-${n%n}; else i=$n; fi
        dir_list_i=$dir_list_i" "$i
    done
    dir_list_i=$(echo -e ${dir_list_i// /\\n} | sort -n)
    smallest=$(echo $dir_list_i | awk '{print $1}')
    largest=$(echo $dir_list_i | awk '{print $NF}')
    step=$(echo "print($(echo $dir_list_i | awk '{print $2}') - ($smallest))" | python)
    echo "Delta from $smallest to $largest with step $step" >> $fname
    echo -e "\nDelta(ratio)\tE(eV)" >> $fname
    
    # still want to get the sorted original dir_list to easily go into the folder named 0.01n.
    unset dir_list
    for i in $dir_list_i
    do
        if [[ "$i" == -* ]]; then n=${i#-}n; else n=$i; fi
        dir_list=$dir_list" "$n
    done
    for n in $dir_list
    do
        if [[ "$n" == *n ]]; then i=-${n%n}; else i=$n; fi
        echo -e "$i\t$(grep sigma $n/OUTCAR | tail -1 | awk '{print $7;}')" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_detector $n
    done
    
    echo $PWD
    echo '' >> $fname
    if [ "$force_converge_flag" ]; then echo -e "!Force doesn't converge during$force_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    if [ "$entropy_converge_flag" ]; then echo -e "!Entropy doesn't converge during$entropy_converge_flag" | tee -a $fname; echo '' >> $fname; fi
    _Display-fit.py $test_type 4 $((4 + data_line_count)) >> $fname

    grep "R-squared is" $fname
    
fi

echo -e "\nMaximal time per run: \c" >> $fname
largest=$(echo $dir_list | awk '{print $NF}')
cat < $(find $largest -mindepth 1 -name "*$largest*") | grep real | awk '{print $2;}' >> $fname
