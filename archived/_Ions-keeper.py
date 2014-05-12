#!/usr/bin/env python
import numpy as np
import sys

scaling_factor_to_be = eval(sys.argv[1])
# open the file that has ions' position part
f = open('POSCAR','r')
file = f.readlines()
f.close

scaling_factor_in_file = eval(file[1])

# convert the raw form of ions' position into a numpy numeral array
ions_position_raw = []
for i in np.arange(8, 16):
    ions_position_raw.append(file[i].split()[0:3])
ions_position = np.array(ions_position_raw, dtype = np.float)

ions_position_new = (ions_position - 0.5) * scaling_factor_in_file / scaling_factor_to_be + 0.5

ions_position_new_raw = np.array(ions_position_new, dtype=np.str).tolist()
for i in range(len(ions_position_new_raw)):
    ions_position_new_raw[i] = ' '.join(ions_position_new_raw[i])+' T T T\n'

#print(ions_position_new_raw)
file[8:16] = ions_position_new_raw
file[1] = str(scaling_factor_to_be) + '\n'
f = open('POSCAR', 'w')
f.writelines(file)
f.close()
