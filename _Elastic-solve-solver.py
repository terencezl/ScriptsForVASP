#!/usr/local/python/2.7.1/bin/python
# Used by Elastic-solve.sh to solve linear equations to get elastic consts
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# _Elastic-solve_solver.py cryst_sys volumn_of_primitive_cell input_data original/alternative

import sys
import numpy as np
import os
import pickle

dirname = os.path.dirname(os.path.realpath(__file__))
cryst_sys = sys.argv[1]
Vpcell = float(sys.argv[2])
econst_input = np.array(eval(sys.argv[3])) * 160.2 / Vpcell

if cryst_sys == 'cubic':
    coeff_matrix = np.array([[1/3., 2/3., 0],
                             [1, -1, 0],
                             [0, 0, 1/2.]])
    if sys.argv[4] == 'alternative':
        f = open(dirname + '/_Elastic-solve-solver-alternative-matrix-cubic.dat', 'rU')
        coeff_matrix = np.array(pickle.load(f))
        f.close()
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC12 is %f\nC44 is %f" % (result[0], result[1], result[2]))

#if cryst_sys == 'cubic_A':
#    coeff_matrix = np.array([[1/2., 0, 0, 0, 0, 0, 0, 0, 0],
#                             [1, 1, 0, 0, 0, 0, 0, 0, 0],
#                             [0, 0, 6, 0, 0, 0, 0, 0, 0],
#                             [0, 0, 0, 1/6., 0, 0, 0, 0, 0],
#                             [0, 0, 0, 1/3., 1, 0, 0, 0, 0],
#                             [0, 0, 0, 1/2., 3, 1, 0, 0, 0],
#                             [0, 0, 0, 1/6., 0, 0, 2, 0, 0],
#                             [0, 0, 0, 1/6., 0, 0, 0, 2, 0],
#                             [0, 0, 0, 0, 0, 0, 0, 0, 8]])
#
#    result = np.linalg.solve(coeff_matrix, econst_input)
#    print("C11 is %f\nC12 is %f\nC44 is %f\
#           \nC111 is %f\nC112 is %f\nC123 is %f\nC144 is %f\nC166 is %f\nC456 is %f"
#               % (result[0], result[1], result[2], result[3], result[4], result[5], result[6], result[7], result[8]))
#    
if cryst_sys == 'tetragonal':
    coeff_matrix = np.array([[1/2., 0, 0, 0, 0, 0],
                             [0, 1/2., 0, 0, 0, 0],
                             [0, 0, 2, 0, 0, 0],
                             [5/2., 1/2., 0, -2, -1, 0],
                             [1, 2, 0, 1, -4, 0],
                             [1, 2, 0, 1, -4, 2]])
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC33 is %f\nC44 is %f\nC12 is %f\nC13 is %f\nC66 is %f"
            % (result[0], result[1], result[2], result[3], result[4], result[5]))

if cryst_sys == 'orthorhombic':
    coeff_matrix = np.array([[1/2., 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 1/2., 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1/2., 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 1/2., 1, 0, 0, 0, 0],
                             [0, 0, 0, 0, 1/2., 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 1/2., 0, 0, 0],
                             [2, 1/2., 1/2., 0, 0, 0, -2, -2, 1],
                             [1/2., 2, 1/2., 0, 0, 0, -2, 1, -2],
                             [1/2., 1/2., 2, 0, 0, 0, 1, -2, -2]])
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC22 is %f\nC33 is %f\nC44 is %f\nC55 is %f\nC66 is %f\nC12 is %f\nC13 is %f\nC23 is %f"
             % (result[0], result[1], result[2], result[3], result[4], 
                result[5], result[6], result[7], result[8]))
