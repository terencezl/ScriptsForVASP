#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

N_steps = 301
Nth_atom = int(sys.argv[1])
#Name_atom = 'Si'
f = open('DOSCAR','rU')
list = []
for line in f:
    list.append(line[0:-1].split())
f.close()

# The projected DOS
Fermi_E = float(list[5][3])
E = []; dos_s = []; dos_p = []; dos_d = []
for n_s in range(0, N_steps):
    E.append(float(list[6+(N_steps+1)*Nth_atom+n_s][0]) - Fermi_E)
    dos_s.append(float(list[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos_p.append(float(list[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos_d.append(float(list[6+(N_steps+1)*Nth_atom+n_s][3]))
plt.plot(E, dos_s, label='s')
plt.plot(E, dos_p, label='p')
plt.plot(E, dos_d, label='d')

plt.legend()
plt.title('Local DOS')
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
#plt.savefig('DOS.png')
plt.show()
