#!/usr/bin/python
# Used by Cellinfo.sh to solve equations

import sys
import numpy as np

task = sys.argv[1]
Vpcell = float(sys.argv[2])

if task =='rwigs':
    N_atoms = eval(sys.argv[3])
    r_input = np.array(eval(sys.argv[4]))
    V_raw = 0
    for i in range(len(r_input)):
        V_raw += 4/3. * np.pi * (N_atoms[i]*r_input[i]**3)
    ratio = (Vpcell/V_raw)**(1/3.)
    r = r_input * ratio
#    print("You should use {} {} as your RWIGS in INCAR to get 100% filling.".format(r[0], r[1]))    # We are dealing with Python 2.7.1 here!
    print("You should use %f %f as your RWIGS in INCAR to get 100%% filling." % (r[0], r[1]))
