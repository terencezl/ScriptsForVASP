#!/usr/local/python/2.7.1/bin/python
import numpy as np
import sys

C_N = a = 1.5
C_H = b = 1.09
N_H = c = 1.01
theta = 20/180. * np.pi
phi = 30/180. * np.pi

C = [0, -a/2, 0]
N = [0, a/2, 0]
H_C_1 = [0, -a/2 - b*np.sin(theta), -b*np.cos(theta)]
H_C_2 = [b*np.cos(theta)*np.cos(phi), -a/2 - b*np.sin(theta), b*np.cos(theta)*np.sin(phi)]
H_C_3 = [-b*np.cos(theta)*np.cos(phi), -a/2 - b*np.sin(theta), b*np.cos(theta)*np.sin(phi)]
H_N_1 = [0, a/2 + c*np.sin(theta), c*np.cos(theta)]
H_N_2 = [c*np.cos(theta)*np.cos(phi), a/2 + c*np.sin(theta), -c*np.cos(theta)*np.sin(phi)]
H_N_3 = [-c*np.cos(theta)*np.cos(phi), a/2 + c*np.sin(theta), -c*np.cos(theta)*np.sin(phi)]

IONS = np.array([C, N, H_C_1, H_C_2, H_C_3, H_N_1, H_N_2, H_N_3])

for i in IONS:
  for j in i:
    print j/10+0.5,
  print('')
