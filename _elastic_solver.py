#!/usr/bin/env python

import sys
import numpy as np

cryst_sys = sys.argv[1]
Vpcell = float(sys.argv[2])
econst_input = np.array(eval(sys.argv[3]))

if cryst_sys == 'cubic':
    coeff_matrix = np.array([[1 / 3., 2 / 3., 0],
                             [1, -1, 0],
                             [0, 0, 1 / 2.]])
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC12 is %f\nC44 is %f" % (result[0], result[1], result[2]))

if cryst_sys == 'tetragonal':
    coeff_matrix = np.array([[1 / 2., 0, 0, 0, 0, 0],
                             [0, 1 / 2., 0, 0, 0, 0],
                             [0, 0, 2, 0, 0, 0],
                             [5 / 2., 1 / 2., 0, -2, -1, 0],
                             [1, 2, 0, 1, -4, 0],
                             [1, 2, 0, 1, -4, 2]])
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC33 is %f\nC44 is %f\nC12 is %f\nC13 is %f\nC66 is %f"
          % (result[0], result[1], result[2], result[3], result[4], result[5]))

if cryst_sys == 'orthorhombic':
    coeff_matrix = np.array([[1 / 2., 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 1 / 2., 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1 / 2., 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 1 / 2., 1, 0, 0, 0, 0],
                             [0, 0, 0, 0, 1 / 2., 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 1 / 2., 0, 0, 0],
                             [2, 1 / 2., 1 / 2., 0, 0, 0, -2, -2, 1],
                             [1 / 2., 2, 1 / 2., 0, 0, 0, -2, 1, -2],
                             [1 / 2., 1 / 2., 2, 0, 0, 0, 1, -2, -2]])
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC22 is %f\nC33 is %f\nC44 is %f\nC55 is %f\nC66 is %f\nC12 is %f\nC13 is %f\nC23 is %f"
          % (result[0], result[1], result[2], result[3], result[4],
             result[5], result[6], result[7], result[8]))
