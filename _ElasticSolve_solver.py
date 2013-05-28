#!/usr/local/python/2.7.1/bin/python
# Used by ElasticSolve.sh to solve linear equations to get elastic consts
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# _ElasticSolve_solver.py Task VolumnOfPrimitiveCell InputData

import sys
import numpy as np

task = sys.argv[1]
Vpcell = float(sys.argv[2])
econst_input = np.array(eval(sys.argv[3])) * 160.2 / Vpcell

if task == 'cubic':
    coeff_matrix = np.array([[3/2., 3, 0],
                             [3, -3, 0],
                             [0, 0, 1/2.]])
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC12 is %f\nC44 is %f" % (result[0], result[1], result[2]))

if task == 'tetragonal':
    coeff_matrix = np.array([[1/2., 0, 0, 0, 0, 0],
                             [0, 1/2., 0, 0, 0, 0],
                             [0, 0, 2, 0, 0, 0],
                             [5/2., 1/2., 0, -2, -1, 0],
                             [1, 2, 0, 1, -4, 0],
                             [1, 2, 0, 1, -4, 2]])
    result = np.linalg.solve(coeff_matrix, econst_input)
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
    result = np.linalg.solve(coeff_matrix, econst_input)
    print("C11 is %f\nC22 is %f\nC33 is %f\nC44 is %f\nC55 is %f\nC66 is %f\nC12 is %f\nC13 is %f\nC23 is %f"
             % (result[0], result[1], result[2], result[3], result[4], 
                result[5], result[6], result[7], result[8]))
