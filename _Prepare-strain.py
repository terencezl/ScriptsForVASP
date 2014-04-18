#!/usr/bin/env python
# Used by Prepare.sh for altering the real space vectors in POSCAR
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118

import sys
import numpy as np

f = open('POSCAR','r')
file = f.readlines()
f.close()

vectors_raw = []
for n in [2,3,4]:
    vectors_raw.append(file[n].split())

vectors = np.array(vectors_raw, dtype=np.float)
delta = float(sys.argv[3])

if sys.argv[2] == 'cubic' or sys.argv[2] == 'cubic_A':
    if sys.argv[1] == "c11+2c12":    # sequence is a little different. Switching 1 and 2
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1 + delta, 0],
                                          [0, 0, 1 + delta]])
    elif sys.argv[1] == "c11-c12":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1 - delta, 0],
                                          [0, 0, 1 + delta**2 / (1-delta**2)]])
    elif sys.argv[1] == "c44":
        transformation_matrix = np.array([[1, delta/2., 0],
                                          [delta/2., 1, 0],
                                          [0, 0, 1 + delta**2 / (4-delta**2)]])
    elif sys.argv[1] == "2c11+2c12+c44":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                         [0, 1 + delta, delta/2],
                                         [0, delta/2, 1]])
    elif sys.argv[1] == "A1":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]])
    elif sys.argv[1] == "A2":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1 + delta, 0],
                                          [0, 0, 1]])
    elif sys.argv[1] == "A3":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1 + delta, 0],
                                          [0, 0, 1 + delta]])
    elif sys.argv[1] == "A4":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1, delta],
                                          [0, delta, 1]])
    elif sys.argv[1] == "A5":
        transformation_matrix = np.array([[1 + delta, delta, 0],
                                          [delta, 1, 0],
                                          [0, 0, 1]])
    elif sys.argv[1] == "A6":
        transformation_matrix = np.array([[1, delta, delta],
                                          [delta, 1, delta],
                                          [delta, delta, 1]])

elif sys.argv[2] == 'tetragonal':
    if sys.argv[1] == "c11":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]])
    elif sys.argv[1] == "c33":
        transformation_matrix = np.array([[1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1 + delta]])
    elif sys.argv[1] == "c44":
        transformation_matrix = np.array([[1, 0, 0],
                                          [0, 1, delta],
                                          [0, delta, 1]])
    elif sys.argv[1] == "5c11-4c12-2c13+c33":
        transformation_matrix = np.array([[1 + 2 * delta, 0, 0],
                                          [0, 1 - delta, 0],
                                          [0, 0, 1 - delta]])
    elif sys.argv[1] == "c11+c12-4c13+2c33":
        transformation_matrix = np.array([[1 - delta, 0, 0],
                                          [0, 1 - delta, 0],
                                          [0, 0, 1 + 2 * delta]])
    elif sys.argv[1] == "c11+c12-4c13+2c33+2c66":
        transformation_matrix = np.array([[1 + delta, delta, 0],
                                          [delta, 1 + delta, 0],
                                          [0, 0, 1 - 2 * delta]])

elif sys.argv[2] == 'orthorhombic':
    if sys.argv[1] == "c11":
        transformation_matrix = np.array([[1 + delta, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]])
    elif sys.argv[1] == "c22":
        transformation_matrix = np.array([[1, 0, 0],
                                          [0, 1 + delta, 0],
                                          [0, 0, 1]])
    elif sys.argv[1] == "c33":
        transformation_matrix = np.array([[1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1 + delta]])
    elif sys.argv[1] == "c44":
        transformation_matrix = np.array([[1, 0, 0],
                                          [0, 1, delta/2],
                                          [0, delta/2, 1]])
    elif sys.argv[1] == "c55":
        transformation_matrix = np.array([[1, 0, delta/2],
                                          [0, 1, 0],
                                          [delta/2, 0, 1]])
    elif sys.argv[1] == "c66":
        transformation_matrix = np.array([[1, delta/2, 0],
                                          [delta/2, 1, 0],
                                          [0, 0, 1]])
    elif sys.argv[1] == "4c11-4c12-4c13+c22+2c23+c33":
        transformation_matrix = np.array([[1 + 2 * delta, 0, 0],
                                          [0, 1 - delta, 0],
                                          [0, 0, 1 - delta]])
    elif sys.argv[1] == "c11-4c12+2c13+4c22-4c23+c33":
        transformation_matrix = np.array([[1 - delta, 0, 0],
                                          [0, 1 + 2 * delta, 0],
                                          [0, 0, 1 - delta]])
    elif sys.argv[1] == "c11+2c12-4c13+c22-4c23+4c33":
        transformation_matrix = np.array([[1 - delta, delta, 0],
                                          [delta, 1 - delta, 0],
                                          [0, 0, 1 + 2 * delta]])

new_vectors = np.array(np.dot(vectors, transformation_matrix), dtype=np.str).tolist()
for n in range(len(new_vectors)):
    new_vectors[n] = ' '.join(new_vectors[n])+'\n'
file[2:5] = new_vectors
f = open('POSCAR','w')
f.writelines(file)
f.close()
