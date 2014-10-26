#!/usr/bin/env bash

### START CONFIG ###
if [[ -z POTENTIAL_DATABASE ]]; then
    POTENTIAL_DATABASE=$HOME/terencelz/local/POTENTIAL_DATABASE
fi
if [[ -z POTENTIAL_TYPE ]]; then
    POTENTIAL_TYPE=PAW-GGA
fi
element_list_file=INPUT_ELEMENT/elements.txt
### END CONFIG ###

while getopts ":e:c:f" opt; do
    case $opt in
    e)
        element_list_file=$OPTARG
        ;;
    c)
        pot_combo=$OPTARG
        ;;
    f)
        is_override=true
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

shift $(($OPTIND-1))
test_type="$1"

if [[ -z "$pot_combo" ]]; then
    echo "-c potential combination must be specified, comma separated, iterating element as X!"
    exit 1
elif [[ $pot_combo == *X* && ! -s $element_list_file ]]; then
    echo "You specified an iterating element X but no element.dat is found!"
    exit 1
elif [[ $pot_combo != *X* ]]; then
    pot_combo_list=$pot_combo
else
    for element in $(cat $element_list_file); do
        pot_combo_item=${pot_combo//X/"$element"}
        pot_combo_list=$pot_combo_list" "$pot_combo_item
    done
fi

pot_combo_list=$(echo $pot_combo_list)
compound_list=${pot_combo_list//,/}


if [[ "$test_type" == prepare ]]; then
    for pot_combo_item in $pot_combo_list
    do
        compound=${pot_combo_item//,/}
        (
        if [[ -d "$compound" ]]; then
            echo -n "$compound/ contains files. "
            if [[ $is_override ]]; then
                echo "Overriding..."
            else
                echo "Escaping..."
                exit 1
            fi
        fi
        mkdir -p "$compound"/INPUT
        cp -r INPUT_ELEMENT/* "$compound"/INPUT
        rm "$compound"/INPUT/*.dat

        if [[ $pot_combo_item != *,* ]]; then
            eval "cat $POTENTIALS_DIR/$POTENTIAL_TYPE/POTCAR_$pot_combo_item > $compound/INPUT/POTCAR"
        else
            eval "cat $POTENTIALS_DIR/$POTENTIAL_TYPE/POTCAR_{$pot_combo_item} > $compound/INPUT/POTCAR"
        fi

        cd "$compound"/INPUT
        sed -i "/SYSTEM/c SYSTEM = $pot_combo_item" INCAR
        sed -i "1c $pot_combo_item" POSCAR
        cd ../..
        )
    done

else
    n=100
    for compound in $compound_list
    do
        (
        if [[ -d "$compound" ]]; then
            echo "Processing $compound/ ..."
            cd "$compound"
        else
            echo "$compound/ does not exist!"
            exit 1
        fi

        i=$compound
        eval "$@"
        )
        n=$(($n+1))
    done
fi
