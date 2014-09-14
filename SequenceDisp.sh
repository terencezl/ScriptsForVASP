#!/usr/bin/env bash

function enter_dir {
    cd "$directory_name"
    if [[ $? != 0 ]]; then
        echo "Cannot find the directory!"
        exit 1
    fi
    echo "$PWD" | tee $fname
}

function prepare_dir_helper {
    # creates global variable dir_list without annoying trailing slashes.
    # be careful about variable name clashes.
    unset dir_list
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
    if [ -n "$force_max" ]; then
        if [ $(echo "$force_max > 0.04" | bc) == 1 ]; then
            force_not_converged_list="$force_not_converged_list\n$dir $force_max"
        fi
    fi
    local entropy=$(grep "entropy T\*S" $dir/OUTCAR | tail -1 | awk '{print $5}')
    entropy=${entropy#-}
    local atom_sum_expr=$(echo $(sed -n '6p' $dir/POSCAR))
    local atom_sum=$((${atom_sum_expr// /+}))
    if [ $(echo "$entropy > 0.001 * $atom_sum" | bc) == 1 ]; then
        entropy_not_converged_list="$entropy_not_converged_list\n$dir $entropy"
    fi
}

function output_force_entropy {
    # takes in force_not_converged_list, entropy_not_converged_list and fname.
    local force_not_converged_list="$1"
    local entropy_not_converged_list="$2"
    local fname="$3"
    if [ "$force_not_converged_list" ]; then
        echo >> $fname
        echo -e "!Force doesn't converge during""$force_not_converged_list" | tee -a $fname
    fi
    if [ "$entropy_not_converged_list" ]; then
        echo >> $fname
        echo -e "!Entropy doesn't converge during""$entropy_not_converged_list" | tee -a $fname
    fi
}

function argparse {
    while getopts ":q:" opt; do
        case $opt in
        q)
            equi_relax=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
      esac
    done
}

if [[ "$1" == */ ]]; then directory_name=${1%/}; else directory_name=$1; fi
test_type="${directory_name%%_*}"
fname="$test_type"_output.txt
data_line_count=0
shift 1

if [[ $test_type == "entest" || $test_type == "kptest" ]]; then
    enter_dir
    prepare_dir_helper
    # get rid of the possible *-1 dirs.
    for dir in $dir_list
    do
        if [[ "$dir" != *-1 ]]; then dir_list_enkp=$dir_list_enkp" "$dir; fi
    done
    dir_list=$(sort_list "$dir_list_enkp")

    prepare_header_helper "$dir_list"
    # echo some headers to file.
    if [ $test_type == "entest" ]; then
        echo -e "ENCUT from $smallest to $largest interval $interval" >> $fname
        echo -e "\n   ENCUT(eV)        E(eV)       dE(eV)" >> $fname
    else
        echo -e "nKP from $smallest to $largest interval $interval" >> $fname
        echo -e "\n         nKP        E(eV)       dE(eV)" >> $fname
    fi
    # echo the data in a sorted way to file.
    for dir in $dir_list
    do
        E=$(grep sigma $dir/OUTCAR | tail -1 | awk '{print $7;}')
        if [ $dir == $smallest ]; then
            E_pre=$E
        fi
        dE=$(echo "($E)-($E_pre)" | bc)
        dE=${dE#-}
        python -c "print '{0:12.6f} {1:12.6f} {2:12.6f}'.format($dir, $E, $dE)" >> $fname
        E_pre=$E
        data_line_count=$(($data_line_count + 1))
    done
    # write the results to file.
    _sequence_fit_plot.py $test_type 4 $data_line_count >> $fname


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

    if [[ $? != 0 ]]; then
        echo "Something is wrong. Terminated."
        exit 1
    fi

    # grep some info to screen.
    grep "!Equilibrium point is out of the considered range!" $fname
    grep "R-squared =" $fname
    grep "Equilibrium scaling constant =" $fname
    grep "B0 =" $fname
    grep "B0' =" $fname
    grep "Total energy =" $fname

    # copy the CONTCAR from the run closest to equilibrium to INPUT/POSCAR.
    # plug equilibrium scaling constant to INPUT/POSCAR
    cd ..
    if [[ -d INPUT ]]; then
        echo "Replacing INPUT/POSCAR with the closest CONTCAR from lctest and updating the scaling constant..."
        scaling_const=$(grep "Equilibrium scaling constant =" \
                    $directory_name/lctest_output.txt | head -1 | awk '{print $5}')
        closest_sc=$(echo $dir_list | awk '{print $1}')
        for dir in $dir_list; do
            diff_dir=$(echo "$dir - $scaling_const" | bc)
            diff_dir=${diff_dir#-}
            diff_sc=$(echo "$closest_sc - $scaling_const" | bc)
            diff_sc=${diff_sc#-}
            if [[ $(echo "$diff_dir < $diff_sc" | bc) == 1 ]]; then closest_sc=$dir; fi
        done
        if [[ -s $directory_name/$closest_sc/CONTCAR ]]; then
            cp $directory_name/$closest_sc/CONTCAR INPUT/POSCAR
            sed -i "2c $scaling_const" INPUT/POSCAR
            cd $directory_name
        else
            echo "Empty CONTCAR. Aborted. Check your runs!"
        fi
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
    equi_relax="equi-relax"
    argparse "$@"
    # some pre-process of the dirs to save computation time.
    enter_dir
    prepare_dir_helper
    for dir in $dir_list; do
        if [[ ! $dir == 0.000 && $dir == *n && ! -d ${dir%n} ]]; then
            cp -rl $dir ${dir%n}
        elif [[ ! $dir == 0.000 && ! $dir == *n && ! -d ${dir}n ]]; then
            cp -rl $dir ${dir}n
        fi
    done
    if [[ ! -d 0.000 ]]; then
        if [[ -d ../$equi_relax ]]; then
            cp -r ../$equi_relax 0.000
        else
            echo "$equi_relax does not exist. It is important because it's used to be 0.000 for the elastic runs."
        fi
    fi

    # get the dir_list and sorted dir_list_minus_sign!
    prepare_dir_helper
    dir_list_minus_sign=$(make_dir_list_sortable "$dir_list")
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
    time_cost=$(grep real $dir/*$dir* 2> /dev/null | awk '{print $2;}')
    if [[ $time_cost ]]; then
        echo $dir' '$time_cost >> $fname
    fi
done
