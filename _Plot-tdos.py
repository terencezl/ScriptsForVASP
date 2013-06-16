#!/usr/local/python/2.7.1/bin/python
import sys
#import numpy as np
import matplotlib.pyplot as plt

N_steps = 301
fname = sys.argv[1]
metal = sys.argv[2]
cryst_struct = sys.argv[3]

f = open(fname,'rU')
list = [];
for line in f:
    list.append(line[0:-1].split())
f.close()

# The total DOS and the integrated DOS
E = []; dos_tot = []; dos_int = []
Fermi_E = float(list[5][3])
for i in range( 6, 6 + N_steps ):
    E.append(float(list[i][0]) - Fermi_E)
    dos_tot.append(float(list[i][1]))
    dos_int.append(float(list[i][2]))
plt.plot(E,dos_tot, label="Total")
#plt.plot(E,dos_int, label="Integrated")

plt.legend()
#plt.axis([0, , , ])
#plt.title("GaAs Total DOS")
plt.xlabel('Energy (eV)')
plt.ylabel('LDOS (State / atom / eV)')
plt.savefig('TDOS.png')
#plt.show()
