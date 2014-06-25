#!/bin/bash
# Prepare.sh entest -s start -e end -i interval -c scaling_const
# Prepare.sh kptest -s start -e end -i interval -c scaling_const
# Prepare.sh lctest -s start -e end -n num_points
# Prepare.sh rttest -s start -e end -n num_points
# Prepare.sh agltest -s start -e end -n num_points -r Ion_rotator_args
# Prepare.sh c11-c12 -y cubic

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

function create_copy_replace {
    mkdir $1 2> /dev/null
    if [[ $? != 0 ]]; then
        if [[ is_override ]]; then
            echo "$1 already exists. Overriding input files."
        else
            echo "$1 already exists. Escaping input files."
        fi
    else
        echo "Creating $1."
    fi
    cd $1
    cp ../../INPUT/INCAR .
    cp ../../INPUT/POSCAR .
    cp ../../INPUT/POTCAR .
    cp ../../INPUT/KPOINTS .
    cp ../../INPUT/qsub.parallel .
    qsub_replacer qsub.parallel
}

function change_dir_name_with_hyphen {
    for dir in $1
    do
        if [[ $dir == -* ]]; then dir=${dir#-}n; fi
        echo $dir
    done
}

function submission_trigger {
    if [[ is_submit ]]; then
        qsub qsub.parallel
    fi
}


directory_name="$1"
echo "Creating test directory $directory_name..."
if [[ -d $directory_name && "$(ls -A $directory_name)" ]]; then
    echo "Warning: Directory contains files or sub-directories. Existing sub-directories will be overriden."
fi
mkdir "$directory_name" 2> /dev/null
cd "$directory_name"
test_type="${directory_name%%_*}"
echo "Preparing $test_type..."
fname="$test_type"_output.txt
shift 1

while getopts ":s:e:n:i:c:y:r:m:o" opt; do
  case $opt in
    s)
      start=$OPTARG
      ;;
    e)
      end=$OPTARG
      ;;
    n)
      num_points=$OPTARG
      ;;
    i)
      interval=$OPTARG
      ;;
    c)
      scaling_const=$OPTARG
      ;;
    y)
      cryst_sys=$OPTARG
      ;;
    r)
      ions_rotator_args=$OPTARG
      ;;
    m)
      is_submit=true
      echo "-m triggered job submission."
      ;;
    o)
      is_override=true
      echo "-o triggered overriding existing subdirectories."
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


if [[ $test_type == "entest" || $test_type == "kptest" ]]; then
    lc1=$scaling_const
    lc2=$(echo "print($lc1+0.1)" | python)
    for ((dir=$start; dir<=$end; dir=$dir+$interval))
    do
        for i in $dir $dir-1
        do
            (
            create_copy_replace $i
            if [ $test_type == "entest" ]; then
                sed -i "s/.*ENCUT.*/ENCUT = $dir/g" INCAR
            else
                sed -i "4c $dir $dir $dir" KPOINTS
            fi
            if [ $i == $dir ]; then
                sed -i "2c $lc1" POSCAR
            else
                sed -i "2c $lc2" POSCAR
            fi
            submission_trigger
            )
        done
    done

elif [[ $test_type == "lctest" ]]; then
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($start,$end,$num_points): print('{0:.3f}'.format(i))" | python)
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        sed -i "2c $dir" POSCAR
        submission_trigger
        )
    done

elif [[ $test_type == "rttest" ]]; then
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($start,$end,$num_points): print('{0:.3f}'.format(i))" | python)
    dir_list=$(change_dir_name_with_hyphen "$dir_list")
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        if [[ "$dir" == *n ]]; then dir=-${dir%n}; fi
        sed -i "s/@R@/$dir/g" POSCAR
        submission_trigger
        )
    done

elif [[ $test_type == "agltest" ]]; then
    dir_list=$(echo -e "import numpy as np\nfor i in np.linspace($start,$end,$num_points): print('{0:.3f}'.format(i))" | python)
    dir_list=$(change_dir_name_with_hyphen "$dir_list")
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        if [[ "$dir" == *n ]]; then dir=-${dir%n}; fi
        Ions_rotator.py "$ions_rotator_args" -a $dir
        submission_trigger
        )
    done

elif [[ $test_type == *c[1-9][1-9]* ]]; then
    if [ $test_type == c44 ]; then
        dir_list="0.020 0.035 0.050n"
    elif [ $test_type == c11-c12 ]; then
        dir_list="0.020 0.030 0.040n"
    fi
    for dir in $dir_list
    do
        (
        create_copy_replace $dir
        if [[ "$dir" == *n ]]; then dir=-${dir%n}; fi
        _prepare_strain.py $test_type $cryst_sys $dir
        submission_trigger
        )
    done

else
    echo "Specify what you are going to test!" >&2
    cd ..
    rm -r "$directory_name"
    exit 1
fi
