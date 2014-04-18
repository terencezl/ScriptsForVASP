#!/bin/bash
# Use: In the top working directory

rename lctest lctest_M lctest || exit 1
scaling_factor_from_M_fit=$(printf %.2f $(grep "Equilibrium scaling factor is" lctest_M/lctest_output.txt | head -1 | awk '{print $5}'))
scaling_factor_start=$(echo "$scaling_factor_from_M_fit - 0.02" | bc)
scaling_factor_end=$(echo "$scaling_factor_from_M_fit + 0.03" | bc)
#Prepare.sh lctest $scaling_factor_start $scaling_factor_start 0.01
#Prepare.sh lctest $scaling_factor_end $scaling_factor_end 0.01
#Fire.sh lctest $scaling_factor_start $scaling_factor_start 0.01
#Fire.sh lctest $scaling_factor_end $scaling_factor_end 0.01
Prep-fire.sh lctest $scaling_factor_start $scaling_factor_end 0.01
#Fire.sh lctest $scaling_factor_start $scaling_factor_end 0.01
