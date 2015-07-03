from __future__ import division
import os
import numpy as np
#import pylab as pl

#Define atom types
Cat = "Ni"; Ani = "S"
LS = "XX" #lattice constant

#Define vectors
a1 = [0, 0.5, 0.5]
a2 = [0.5, 0, 0.5]
a3 = [0.5, 0.5, 0]

#Define internal parameter for spinel structure
pr1 = -0.135

#for 8a (0,0,0)

x = [float(0)]*14;y = [float(0)]*14; z = [float(0)]*14; u = [float(0)]*14; v = [float(0)]*14; w = [float(0)]*14
np.array(x);np.array(y);np.array(z);np.array(u);np.array(v);np.array(w)
x[0] = 0.; y[0] = 0.; z[0] = 0.
x[1] = 0.75; y[1] = 0.25; z[1] = 0.75
x[2] = 0.625; y[2] = 0.625; z[2] = 0.625
x[3] = 0.375; y[3] = 0.875; z[3] = 0.125
x[4] = 0.875; y[4] = 0.125; z[4] = 0.375
x[5] = 0.125; y[5] = 0.375; z[5] = 0.875
x[6] = pr1; y[6] = pr1; z[6] = pr1
x[7] = 0.5 - pr1; y[7] = 0.5 + pr1; z[7] = -pr1
x[8] = 0.75 + pr1; y[8] = 0.25 + pr1; z[8] = 0.75 - pr1
x[9] = 0.25 + pr1; y[9] = 0.75 - pr1; z[9] = 0.75 + pr1
x[10] = -pr1; y[10] = 0.5 - pr1; z[10] = 0.5 + pr1
x[11] = 0.5 + pr1; y[11] = -pr1; z[11] = 0.5 - pr1
x[12] = 0.25 - pr1; y[12] = 0.25 - pr1; z[12] = 0.25 - pr1
x[13] = 0.75 - pr1; y[13] = 0.75 + pr1; z[13] = 0.25 + pr1

i=0
while (i<14):
    u[i] = y[i] + z[i] - x[i]
    v[i] = x[i] + z[i] - y[i]
    w[i] = x[i] + y[i] - z[i]
    if u[i]<0:
        u[i] = u[i] + 1
    if u[i]>1:
        u[i] = u[i] - 1
    if v[i]<0:
        v[i] = v[i] + 1
    if v[i]>1:
        v[i] = v[i] - 1
    if w[i]<0:
        w[i] = w[i] + 1
    if w[i]>1:
        w[i] = w[i] - 1     
    i+=1

#np.savetxt("POSCAR_spinel", np.column_stack((u, v, w)))

#write POSCAR file
pos = open('POSCAR', 'wb')
s = str(Cat) + "3 " + str(Ani) + "4 Spinel"
pos.write(s)
pos.write('\n')
pos.write(LS)
pos.write('\n')

s = str.format("{0:.5f}", a1[0]) + "    " + str.format("{0:.5f}", a1[1]) + "    " + str.format("{0:.5f}", a1[2])
pos.write(s)
pos.write('\n')
s = str.format("{0:.5f}", a2[0]) + "    " + str.format("{0:.5f}", a2[1]) + "    " + str.format("{0:.5f}", a2[2])
pos.write(s);
pos.write('\n')
s = str.format("{0:.5f}", a3[0]) + "    " + str.format("{0:.5f}", a3[1]) + "    " + str.format("{0:.5f}", a3[2])
pos.write(s)
pos.write('\n')
s = str(Cat) + "    " +str(Ani)
pos.write(s)
pos.write('\n')
s = "6    8"
pos.write(s)
pos.write('\n')
s = "Direct"
pos.write(s)
pos.write('\n')
i = 0
while (i<14):
    s = str.format("{0:.11f}", u[i]) + "  " + str.format("{0:.11f}", v[i]) + "  " + str.format("{0:.11f}", w[i])
    pos.write(s)
    pos.write('\n')
    i+= 1
pos.close
