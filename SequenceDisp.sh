#!/usr/bin/env bash

function enter_dir {
    cd "$directory_name"
    if [[ $? != 0 ]]; then
        echo "Cannot find the directory!"
        exit 1
    fi
}

function prepare_dir_helper {
    # creates global variable dir_list without annoying trailing slashes.
    # be careful about variable name clashes.
    echo "$PWD" | tee $fname
    for dir in $(ls -F)
    do
        if [[ "$dir" == */ ]]; then dir_list=$dir_list" "${dir%/}; fi
    done
}

function make_dir_list_sortable {
    # takes one dir_list.
    # echo returns the sortable list, with *n turned to -*.
    local dir_list="$1"
    local dir_list_minus_sign
    for dir in $dir_list
    do
        if [[ "$dir" == *n ]]; then dir_minus_sign=-${dir%n}; else dir_minus_sign=$dir; fi
        dir_list_minus_sign=$dir_list_minus_sign" "$dir_minus_sign
    done
    echo $dir_list_minus_sign
}

function sort_list {
    # takes one sortable dir_list.
    # echo returns the sorted dir_list.
    local dir_list_sortable="$1"
    local dir_list_sorted
    dir_list_sortable="${dir_list_sortable// /\\n}"
    dir_list_sorted=$(echo -e "$dir_list_sortable" | sort -n)
    echo $dir_list_sorted
}

function prepare_header_helper {
    # creates global variables smallest, largest, second_small and interval.
    # be careful about variable name clashes.
    local dir_list_sorted="$1"
    smallest=$(echo $dir_list_sorted | awk '{print $1}')
    largest=$(echo $dir_list_sorted | awk '{print $NF}')
    second_small=$(echo $dir_list_sorted | awk '{print $2}')
    interval=$(python -c "print('{0:.3f}'.format($second_small - $smallest))")
}

function force_entropy_not_converged_detecting_helper {
    # takes one dir
    # creates global variables force_not_converged_list and entropy_not_converged_list
    # be careful about name clashes.
    local dir="$1"
    local force_max=$(grep "FORCES:" $dir/OUTCAR | tail -1 | awk '{print $5}')
    force_max=${force_max#-}
    if [ $(echo "$force_max > 0.04" | bc) == 1 ]; then
        force_not_converged_list="$force_not_converged_list$1 $force_max\n"
    fi

    local entropy=$(grep "entropy T\*S" $dir/OUTCAR | tail -1 | awk '{print $5}')
    entropy=${entropy#-}
    local atom_sum_expr=$(echo $(sed -n '6p' $dir/POSCAR))
    local atom_sum=$((${atom_sum_expr// /+}))
    if [ $(echo "$entropy > 0.001 * $atom_sum" | bc) == 1 ]; then
        entropy_not_converged_list="$entropy_not_converged_list$1 $entropy\n"
    fi
}

function output_force_entropy {
    # takes in force_not_converged_list, entropy_not_converged_list and fname.
    local force_not_converged_list="$1"
    local entropy_not_converged_list="$2"
    local fname="$3"
    if [ "$force_not_converged_list" ]; then
        echo -e "!Force doesn't converge during\n$force_not_converged_list" | tee -a $fname
    fi
    if [ "$entropy_not_converged_list" ]; then
        echo -e "!Entropy doesn't converge during\n$entropy_not_converged_list" | tee -a $fname
    fi
}


if [[ "$1" == */ ]]; then directory_name=${1%/}; else directory_name=$1; fi
test_type="${directory_name%%_*}"
fname="$test_type"_output.txt
data_line_count=0


if [[ $test_type == "entest" || $test_type == "kptest" ]]; then
    enter_dir
    prepare_dir_helper
    # get rid of the *-1 dirs.
    for dir in $dir_list
    do
        if [[ "$dir" != *-1 ]]; then dir_list_enkp=$dir_list_enkp" "$dir; fi
    done
    dir_list=$(sort_list "$dir_list_enkp")
    prepare_header_helper "$dir_list"
    lc1=$(sed -n '2p' $smallest/POSCAR)
    lc2=$(echo "$lc1+0.1" | bc)
    # echo some headers to file.
    if [ $test_type == "entest" ]; then
        echo -e "ENCUT from $smallest to $largest interval $interval" >> $fname
        echo -e "LC1 = $lc1, LC2 = $lc2 (directories with '-1')" >> $fname
        echo -e "\n   ENCUT(eV)    E_LC1(eV)   dE_LC1(eV)    E_LC2(eV)   dE_LC2(eV) DE=E_LC2-E_LC1(eV)      dDE(eV)" >> $fname
    else
        echo -e "nKP from $smallest to $largest interval $interval" >> $fname
        echo -e "LC1 = $lc1, LC2 = $lc2 (directories with '-1')" >> $fname
        echo -e "\n         nKP    E_LC1(eV)   dE_LC1(eV)    E_LC2(eV)   dE_LC2(eV) DE=E_LC2-E_LC1(eV)      dDE(eV)" >> $fname
    fi
    # echo the data in a sorted way to file.
    for dir in $dir_list
    do
        E_LC1=$(grep sigma $dir/OUTCAR | tail -1 | awk '{print $7;}')
        E_LC2=$(grep sigma $dir-1/OUTCAR | tail -1 | awk '{print $7;}')
        DE=$(echo "($E_LC2)-($E_LC1)" | bc)
        if [ $dir == $smallest ]; then
            DE_pre=$DE
            E_LC1_pre=$E_LC1
            E_LC2_pre=$E_LC2
        fi
        dDE=$(echo "($DE)-($DE_pre)" | bc)
        dE_LC1=$(echo "($E_LC1)-($E_LC1_pre)" | bc)
        dE_LC2=$(echo "($E_LC2)-($E_LC2_pre)" | bc)
        dDE=${dDE#-}
        dE_LC1=${dE_LC1#-}
        dE_LC2=${dE_LC2#-}
        python -c "print '{0:12.6f} {1:12.6f} {2:12.6f} {3:12.6f} {4:12.6f} {5:18.6f} {6:12.6f}'.format(\
                                    $dir, $E_LC1, $dE_LC1, $E_LC2, $dE_LC2, $DE, $dDE)" >> $fname
        # set DE for this cycle as DE_pre for the next cycle to get the next dDE
        DE_pre=$DE
        E_LC1_pre=$E_LC1
        E_LC2_pre=$E_LC2
        data_line_count=$(($data_line_count + 1))
    done
    # write the results to file.
    _sequence_fit_plot.py $test_type 5 $data_line_count >> $fname


elif [[ $test_type == "lctest" ]]; then
    enter_dir
    prepare_dir_helper
    dir_list=$(sort_list "$dir_list")
    prepare_header_helper "$dir_list"
    # echo some headers to file.
    echo -e "Scaling constant from $smallest to $largest interval $interval" >> $fname
    echo -e "\nScalingConst(Ang) Volume(Ang^3)         E(eV)" >> $fname
    # echo the data in a sorted way to file.
    for dir in $dir_list
    do
        Vpcell=$(cat $dir/OUTCAR | grep 'volume of cell' | tail -1 | awk '{print $5;}')
        energy=$(grep sigma $dir/OUTCAR | tail -1 | awk '{print $7;}')
        python -c "print '{0:17.6f} {1:13.6f} {2:13.6f}'.format($dir, $Vpcell, $energy)" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_not_converged_detecting_helper $dir
    done
    # echo runs with not converged force or entropy to file and screen.
    output_force_entropy "$force_not_converged_list" "$entropy_not_converged_list" "$fname"
    # fit the data and write the results to file.
    _sequence_fit_plot.py $test_type 4 $data_line_count>> $fname
    # grep some info to screen.
    grep "!Equilibrium point is out of the considered range!" $fname
    grep "R-squared =" $fname
    grep "Equilibrium scaling constant =" $fname
    grep "B0 =" $fname
    grep "B0' =" $fname
    grep "Total energy =" $fname

    # copy the CONTCAR from the run closest to equlibrium to INPUT/POSCAR.
    # plug equlibirum scaling constant to INPUT/POSCAR
    if [[ -d INPUT ]]; then
        cd ..
        scaling_const=$(grep "Equilibrium scaling constant =" \
                    $directory_name/lctest_output.txt | head -1 | awk '{print $5}')
        closest_lc=$(echo $dir_list | awk '{print $1}')
        for dir in $dir_list; do
            diff_dir=$(echo "$dir - $scaling_const" | bc)
            diff_dir=${diff_dir#-}
            diff_lc=$(echo "$closest_lc - $scaling_const" | bc)
            diff_lc=${diff_lc#-}
            if [[ $(echo "$diff_dir < $diff_lc" | bc) == 1 ]]; then closest_lc=$dir; fi
        done
        cp $directory_name/$closest_lc/CONTCAR INPUT/POSCAR
        sed -i "2c $scaling_const" INPUT/POSCAR
        cd $directory_name
    fi


elif [[ $test_type == "rttest" ]]; then
    enter_dir
    prepare_dir_helper
    dir_list=$(sort_list "$dir_list")
    prepare_header_helper "$dir_list"
    # echo some headers to file.
    echo -e "Ratio from $smallest to $largest interval $interval" >> $fname
    echo -e "\n        Ratio        Volume         E(eV)" >> $fname
    # echo the data in a sorted way to file.
    for dir in $dir_list
    do
        Vpcell=$(cat $dir/OUTCAR | grep 'volume of cell' | tail -1 | awk '{print $5;}')
        energy=$(grep sigma $dir/OUTCAR | tail -1 | awk '{print $7;}')
        python -c "print '{0:13.6f} {1:13.6f} {2:13.6f}'.format($dir, $Vpcell, $energy)" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_not_converged_detecting_helper $dir
    done
    # echo runs with not converged force or entropy to file and screen.
    output_force_entropy "$force_not_converged_list" "$entropy_not_converged_list" "$fname"
    # fit the data and write the results to file.
    _sequence_fit_plot.py $test_type 4 $data_line_count>> $fname
    # grep some info to screen.
    grep "!Equilibrium point = out of the considered range!" $fname
    grep "R-squared =" $fname
    grep "Equilibrium ratio =" $fname
    grep "B0 =" $fname
    grep "B0' =" $fname
    grep "Total energy =" $fname


elif [[ $test_type == "agltest" ]]; then
    enter_dir
    prepare_dir_helper
    # get the dir_list and sorted dir_list_minus_sign!
    dir_list_minus_sign=$(make_dir_list_sortable $dir_list)
    dir_list_minus_sign=$(sort_list "$dir_list_minus_sign")
    prepare_header_helper "$dir_list_minus_sign"
    # echo some headers to file.
    echo -e "Angle from $smallest to $largest interval $interval" >> $fname
    echo -e "\n        Angle         E(eV)" >> $fname
    # echo the data in a sorted way to file.
    for dir_minus_sign in $dir_list_minus_sign
    do
        if [[ "$dir_minus_sign" == -* ]]; then dir=${dir_minus_sign#-}n; else dir=$dir_minus_sign; fi
        energy=$(grep sigma $dir/OUTCAR | tail -1 | awk '{print $7;}')
        python -c "print '{0:13.6f} {1:13.6f}'.format($dir_minus_sign, $energy)" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_not_converged_detecting_helper $dir
    done
    # echo runs with not converged force or entropy to file and screen.
    output_force_entropy "$force_not_converged_list" "$entropy_not_converged_list" "$fname"
    # write the results to file.
    _sequence_fit_plot.py $test_type 4 $data_line_count>> $fname


elif [[ $test_type == "equi-relax" ]]; then
    cd "$directory_name"
    if [[ $? != 0 ]]; then
        echo "Cannot find the directory!"
        exit 1
    fi
    cp CONTCAR ../INPUT/POSCAR
    exit 0


elif [[ $test_type == *c[1-9][1-9]* ]]; then
    # some pre-process of the dirs to save computation time.
    enter_dir
    function duplicate {
        dir_from=$1
        if [[ "$dir_from" == *n ]]; then dir_to=${dir_from%n}; else dir_to="$dir_from"n; fi
        if [ -d $dir_to ]; then rm -r $dir_to; fi
        cp -rl $dir_from $dir_to
    }

    if [[ $test_type == c44 ]]; then
        duplicate 0.020
        duplicate 0.035
        duplicate 0.050n
    elif [[ $test_type == c11-c12 ]]; then
        duplicate 0.020
        duplicate 0.030
        duplicate 0.040n
    fi

    if [ -d 0.000 ]; then rm -r 0.000; fi
    cp -r ../../equi-relax 0.000

    # get the dir_list and sorted dir_list_minus_sign!
    prepare_dir_helper
    dir_list_minus_sign=$(make_dir_list_sortable $dir_list)
    dir_list_minus_sign=$(sort_list "$dir_list_minus_sign")
    prepare_header_helper "$dir_list_minus_sign"
    # echo some headers to file.
    echo "Delta from $smallest to $largest with interval $interval" >> $fname
    echo -e "\nDelta(ratio)         E(eV)" >> $fname
    # echo the data in a sorted way to file. 
    for dir_minus_sign in $dir_list_minus_sign
    do
        if [[ "$dir_minus_sign" == -* ]]; then dir=${dir_minus_sign#-}n; else dir=$dir_minus_sign; fi
        energy=$(grep sigma $dir/OUTCAR | tail -1 | awk '{print $7;}')
        python -c "print '{0:12.6f} {1:13.6f}'.format($dir_minus_sign, $energy)" >> $fname
        data_line_count=$(($data_line_count + 1))
        force_entropy_not_converged_detecting_helper $dir
    done
    # echo runs with not converged force or entropy to file and screen.
    output_force_entropy "$force_not_converged_list" "$entropy_not_converged_list" "$fname"
    # fit the data and write the results to file.
    _sequence_fit_plot.py $test_type 4 $data_line_count >> $fname
    # grep some info to screen.
    grep "R-squared =" $fname

else
    echo "Specify what you are going to display!" >&2
    exit 1
fi

echo -e "\ntime cost per run:" >> $fname
for dir in $dir_list; do
    echo -n $dir' ' >> $fname
    grep real $dir/*$dir* | awk '{print $2;}' >> $fname
done
