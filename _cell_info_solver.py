#!/usr/bin/env python
# Used by CellInfo.sh to solve equations
# _cell_info_solver.py test_type volumn_of_primitive_cell input_data(in the form of a python list)

import sys
import numpy as np

test_type = sys.argv[1]
Vpcell = float(sys.argv[2])

if test_type == 'rwigs':
    N_atoms = eval(sys.argv[3])
    r_input = np.array(eval(sys.argv[4]))
    V_raw = 0
    for i in range(len(r_input)):
        V_raw += 4 / 3. * np.pi * (N_atoms[i] * r_input[i] ** 3)
    ratio = (Vpcell / V_raw) ** (1 / 3.)
    r = r_input * ratio
    print "You should use ",
    if len(r) == 1:
        sys.stdout.write("%.5f" % r[0])
    else:
        for i in range(len(r)):
            sys.stdout.write("%.5f," % r[i])
    print " as your RWIGS in INCAR to get 100% filling."