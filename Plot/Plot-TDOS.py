#!/usr/bin/python
import sys
#import numpy as np
import matplotlib.pyplot as plt

N_steps = 301
#if len(sys.argv) != 2:
#    N_steps = int(sys.argv[-1])

f = open('DOSCAR','rU')
list = [];
for line in f:
    list.append(line[0:-1].split())
f.close()

# The total DOS and the integrated DOS.
E = []; dos_tot = []; dos_int = []
Fermi_E = float(list[5][3])
for i in range( 6, 6 + N_steps ):
    E.append(float(list[i][0]) - Fermi_E)
    dos_tot.append(float(list[i][1]))
    dos_int.append(float(list[i][2]))
plt.plot(E,dos_tot, label="Total")
#plt.plot(E,dos_int, label="Integrated")

#plt.axis([0, , , ])
plt.legend()
plt.title("GaAs Total DOS")
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
#plt.savefig('TDOS.png')
plt.show()
