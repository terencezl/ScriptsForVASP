#!/usr/bin/env python
# Usage: Ions_rotator.py '[ions line number as a list]' '[rotation center vector in direct coords as a list]'
# '[axis directional vector in cartesian coords as a list(no need to be normalized)]'
# angle(counter-clockwise) POSCAR_file_name
# e.g.: rotator.py '[23,24,25,26,27,28]' '[0,0.5,0.75]' '[np.sqrt(2),0,1]' 60 POSCAR

import numpy as np
import sys

np.set_printoptions(suppress=True)

ions_line_number = np.array(eval(sys.argv[1])) - 1  # machine counts from 0
rotation_center_direct_vector = eval(sys.argv[2])

# define the rotation matrix A. Here we use the active rotation scheme, i.e. rotating the object
# axis = [u_x, u_y, u_z] = [np.sqrt(2), 0, 1]/np.sqrt(3)
axis = eval(sys.argv[3])
[u_x, u_y, u_z] = axis / np.linalg.norm(axis)
angle = eval(sys.argv[4]) / 180. * np.pi
A = np.eye(3) * np.cos(angle) + np.array([[0, -u_z, u_y], [u_z, 0, -u_x], [-u_y, u_x, 0]]) \
                                * np.sin(angle) + np.array([[u_x ** 2, u_x * u_y, u_x * u_z],
                                                            [u_x * u_y, u_y ** 2, u_y * u_z],
                                                            [u_x * u_z, u_y * u_z, u_z ** 2]]) * (1 - np.cos(angle))

# open the POSCAR that has ions' position part
with open(sys.argv[5], 'r') as f:
    POSCAR = f.readlines()

# achieve the basis vectors from POSCAR
basis_vectors_raw = []
for i in range(2, 5):
    basis_vectors_raw.append(POSCAR[i].split())
basis_vectors = np.array(basis_vectors_raw, dtype=np.float)

# convert the raw form of ions' position into a numpy numeral array and subtracting the translational vector
ions_position_raw = []
for i in ions_line_number:
    ions_position_raw.append(POSCAR[i].split()[0:3])
ions_position = np.array(ions_position_raw, dtype=np.float) - rotation_center_direct_vector

# transform the coordinates from direct to cartesian
ions_position = np.dot(ions_position, basis_vectors)

# the matrix multiplication
ions_position_new = np.dot(ions_position, A.transpose())

# transform the coordinates from cartesian back to direct
ions_position_new = np.dot(ions_position_new, np.linalg.inv(basis_vectors))

# ions_position_new = ions_position_new / np.cos(angle)

# add the translational vector
ions_position_new = ions_position_new + rotation_center_direct_vector

# display the rotated ion positions
for i in ions_position_new:
    print('{0:.8f} {1:.8f} {2:.8f}'.format(i[0], i[1], i[2]))

ions_position_new_raw = np.array(ions_position_new, dtype=np.str).tolist()
for i in range(len(ions_position_new_raw)):
    ions_position_new_raw[i] = ' '.join(ions_position_new_raw[i]) + '\n'

count = 0
for i in ions_line_number:
    POSCAR[i] = ions_position_new_raw[count]
    count += 1
#f = open('POSCAR', 'w')
with open(sys.argv[5], 'w') as f:
    f.writelines(POSCAR)