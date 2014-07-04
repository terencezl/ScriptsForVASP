#!/usr/bin/env bash
# Use: In the TMN working directory where there is a list of folders by the name of elements

POT_TYPE=PAW-GGA
POTENTIALS_DIR=$HOME/terencelz/local/potential-database

element_list_file=INPUT_ELEMENT/element.dat

while getopts ":e:c:" opt; do
    case $opt in
    e)
        element_list_file=$OPTARG
        ;;
    c)
        pot_combo=$OPTARG
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
        mkdir -p "$compound"/INPUT
        cp -r INPUT_ELEMENT/* "$compound"/INPUT
        rm "$compound"/INPUT/*.dat

        if [[ $pot_combo_item != *,* ]]; then
            eval "cat $POTENTIALS_DIR/$POT_TYPE/POTCAR_$pot_combo_item > $compound/INPUT/POTCAR"
        else
            eval "cat $POTENTIALS_DIR/$POT_TYPE/POTCAR_{$pot_combo_item} > $compound/INPUT/POTCAR"
        fi

        cd "$compound"/INPUT
        sed -i "/SYSTEM/c SYSTEM = $pot_combo_item" INCAR
        sed -i "1c $pot_combo_item" POSCAR
        cd ../..
    done

else
#    test_script=$test_type
#    shift 1
#    if ! type $test_script >/dev/null 2>&1; then
#        echo "Command $test_script does not exist!"
#        exit 1
#    fi
    for compound in $compound_list
    do
        (
#        if ! cd "$compound"; then
#            echo "$compound directory does not exist!"
#            exit 1
#        fi
        cd "$compound"
        $@
        )
    done
fi
