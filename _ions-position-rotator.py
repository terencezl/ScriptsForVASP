#!/usr/local/python/2.7.1/bin/python
import numpy as np
import sys
np.set_printoptions(suppress=True)

angle = eval(sys.argv[1]) / 180. * np.pi
#angle = np.arctan(1/np.sqrt(2))

# define the rotation matrix A. Here we use the active rotation scheme, i.e. rotating the object
#axis = [u_x, u_y, u_z] = [-1/np.sqrt(2), 1/np.sqrt(2), 0]
axis = [u_x, u_y, u_z] = [-1, 1, 0]/np.sqrt(2)

#A= np.array([[np.cos(angle)+u_x**2*(1-np.cos(angle)), u_x*u_y*(1-np.cos(angle))-u_z*np.sin(angle), u_x*u_z*(1-np.cos(angle))+u_y*np.sin(angle)],
#             [u_y*u_x*(1-np.cos(angle))+u_z*np.sin(angle), np.cos(angle)+u_y**2*(1-np.cos(angle)), u_y*u_z*(1-np.cos(angle))-u_x*np.sin(angle)],
#             [u_z*u_x*(1-np.cos(angle))-u_y*np.sin(angle), u_z*u_y*(1-np.cos(angle))+u_x*np.sin(angle), np.cos(angle)+u_z**2*(1-np.cos(angle))]])
# an alternative form
A = np.eye(3) * np.cos(angle) + np.array([[0, -u_z, u_y], [u_z, 0, -u_x], [-u_y, u_x, 0]]) * np.sin(angle) + np.array([[u_x**2, u_x*u_y, u_x*u_z], [u_x*u_y, u_y**2, u_y*u_z], [u_x*u_z, u_y*u_z, u_z**2]]) * (1 - np.cos(angle))

# open the file that has ions' position part
f = open('POSCAR','r')
file = f.readlines()
f.close

# convert the raw form of ions' position into a numpy numeral array
ions_position_raw = []
for i in np.arange(8, 16):
    ions_position_raw.append(file[i].split()[0:3])
ions_position = np.array(ions_position_raw, dtype = np.float) - 0.5

# the matrix multiplication
ions_position_new = np.dot(ions_position, A.transpose()) + 0.5

# display the rotated ion positions
#for i in ions_position_prime:
#  for j in i:
#    print j,
#  print('')

ions_position_new_raw = np.array(ions_position_new, dtype=np.str).tolist()
for i in range(len(ions_position_new_raw)):
    ions_position_new_raw[i] = ' '.join(ions_position_new_raw[i])+' T T T\n'

file[8:16] = ions_position_new_raw
f = open('POSCAR', 'w')
f.writelines(file)
f.close()
