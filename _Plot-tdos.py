#!/usr/local/python/2.7.1/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import stineman_interp
#import pdb

with open('DOSCAR','r') as f:
    list = f.readlines()
for i in range(0, len(list)):
    list[i] = list[i].split()

N_steps = int(list[5][2])
Fermi_E = float(list[5][3])

# The total DOS and the integrated DOS
E = []; dos_tot = []; dos_int = []
Fermi_E = float(list[5][3])
for i in range(6, 6 + N_steps):
    E.append(float(list[i][0]) - Fermi_E)
    dos_tot.append(float(list[i][1]))
    dos_int.append(float(list[i][2]))
plt.plot(E, dos_tot)
#plt.plot(E, dos_int, label="Integrated")

table = np.column_stack((E, dos_tot, dos_int))
np.savetxt('TDOS.txt', table, '%.6f', '\t')

slice = dos_int[np.abs(np.array(E) - 0.2).argmin()] - dos_int[np.abs(np.array(E) + 0.2).argmin()]
#slice = dos_tot[np.abs(np.array(E)).argmin()]
np.savetxt('TDOS@Ef.txt', [slice], '%.6f')

#plt.legend()
#plt.axis([-20, 5, 0, 20])
plt.xlabel('Energy (eV)')
plt.ylabel('TDOS (States / Unit Cell / eV)')
plt.savefig('TDOS.png')
#plt.show()
