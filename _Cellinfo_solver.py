#!/usr/bin/python
# Used by Cellinfo.sh to solve equations
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# _Cellinfo_solver.py Task InputData VolumnOfPrimitiveCell

import sys
import numpy as np

task = sys.argv[1]
Vpcell = float(sys.argv[3])
econst_raw = np.array(eval(sys.argv[2])) * 160.2 / Vpcell

if task =='rwigs':
    N = eval(sys.argv[2])
    r_raw = np.array(eval(sys.argv[4]))
    V_raw = 0
    for i in range(len(r_raw)):
        V_raw += 4/3. * np.pi * (N[i]*r_raw[i]**3)
    ratio = (Vpcell/V_raw)**(1/3.)
    r = r_raw * ratio
#    print("You should use {} {} as your RWIGS in INCAR to get 100% filling.".format(r[0], r[1]))    # We are dealing with Python 2.7.1 here!
    print("You should use %f %f as your RWIGS in INCAR to get 100%% filling." % (r[0], r[1]))

if task == 'cubic':
    coeff_matrix = np.array([[3/2., 3, 0],
                             [3, -3, 0],
                             [0, 0, 1/2.]])
    result = np.linalg.solve(coeff_matrix, econst_raw)
#    print("C11 is {}\nC12 is {}\nC44 is {}".format(result[0], result[1], result[2]))
    print("C11 is %f\nC12 is %f\nC44 is %f" % (result[0], result[1], result[2]))

if task == 'tetragonal':
    coeff_matrix = np.array([[1/2., 0, 0, 0, 0, 0],
                             [0, 1/2., 0, 0, 0, 0],
                             [0, 0, 2, 0, 0, 0],
                             [5/2., 1/2., 0, -2, -1, 0],
                             [1, 2, 0, 1, -4, 0],
                             [1, 2, 0, 1, -4, 2]])
    result = np.linalg.solve(coeff_matrix, econst_raw)
#    print("C11 is {}\nC33 is {}\nC44 is {}\nC12 is {}\nC13 is {}\nC66 is {}"
#            .format(result[0], result[1], result[2], result[3], result[4], result[5]))
    print("C11 is %f\nC33 is %f\nC44 is %f\nC12 is %f\nC13 is %f\nC66 is %f"
            % (result[0], result[1], result[2], result[3], result[4], result[5]))

if task == 'orthorhombic':
    coeff_matrix = np.array([[1/2., 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 1/2., 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1/2., 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 1/2., 1, 0, 0, 0, 0],
                             [0, 0, 0, 0, 1/2., 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 1/2., 0, 0, 0],
                             [2, 1/2., 1/2., 0, 0, 0, -2, -2, 1],
                             [1/2., 2, 1/2., 0, 0, 0, -2, 1, -2],
                             [1/2., 1/2., 2, 0, 0, 0, 1, -2, -2]])
    result = np.linalg.solve(coeff_matrix, econst_raw)
#    print("C11 is {}\nC22 is {}\nC33 is {}\nC44 is {}\nC55 is {}\nC66 is {}\nC12 is {}\nC13 is {}\nC23 is {}"
#            .format(result[0], result[1], result[2], result[3], result[4], 
#                result[5], result[6], result[7], result[8]))
    print("C11 is %f\nC22 is %f\nC33 is %f\nC44 is %f\nC55 is %f\nC66 is %f\nC12 is %f\nC13 is %f\nC23 is %f"
             % (result[0], result[1], result[2], result[3], result[4], 
                result[5], result[6], result[7], result[8]))
