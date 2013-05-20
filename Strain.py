#!/usr/bin/python
import sys
import numpy as np
#import matplotlib.pyplot as plt
import pdb
f = open('POSCAR','r')
#file = []
#for line in f:
#    file.append(line.split())
file = f.readlines()
f.close()

vectors_raw = []
for n in [2,3,4]:
    vectors_raw.append(file[n].split())

vectors = np.array(vectors_raw, dtype=np.float)
#pdb.set_trace()
delta = 0.01
delta = float(sys.argv[2])

if sys.argv[1] == "c11-c12":
    transformation_matrix = np.array([[1 + delta, 0, 0], [0, 1 + delta, 0], [0, 0, (1 + delta)**(-2)]])
elif sys.argv[1] == "c11+2c12":
    transformation_matrix = np.array([[1 + delta, 0, 0], [0, 1 + delta, 0], [0, 0, 1 + delta]])
elif sys.argv[1] == "c44":
    transformation_matrix = np.array([[1, delta/2, 0], [delta/2, 1, 0], [0, 0, 1 + delta**2 / (4-delta**2)]])

new_vectors = np.array(np.dot(vectors, transformation_matrix), dtype=np.str).tolist()
for n in range(len(new_vectors)):
    new_vectors[n] = ' '.join(new_vectors[n])+'\n'
file[2:5] = new_vectors
f = open('POSCAR','w')
f.writelines(file)
f.close()


