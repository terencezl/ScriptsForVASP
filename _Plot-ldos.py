#!/usr/local/python/2.7.1/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

with open('DOSCAR','r') as f:
    list = f.readlines()
for i in range(0, len(list)):
    list[i] = list[i].split()

N_steps = int(list[5][2])
Fermi_E = float(list[5][3])

# The projected DOS
Nth_atom = 1
E = []; dos1_s = []; dos1_p = []; dos1_d = []
for n_s in range(0, N_steps):
    E.append(float(list[6+(N_steps+1)*Nth_atom+n_s][0]) - Fermi_E)
    dos1_s.append(float(list[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos1_p.append(float(list[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos1_d.append(float(list[6+(N_steps+1)*Nth_atom+n_s][3]))
plt.plot(E, dos1_s, label= str(Nth_atom) + '_s')
plt.plot(E, dos1_p, label= str(Nth_atom) + '_p')
plt.plot(E, dos1_d, label= str(Nth_atom) + '_d')

Nth_atom = 2
dos2_s = []; dos2_p = []; dos2_d = []
for n_s in range(0, N_steps):
    dos2_s.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos2_p.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos2_d.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][3]))
plt.plot(E, dos2_s, label= str(Nth_atom) + '_s')
plt.plot(E, dos2_p, label= str(Nth_atom) + '_p')
plt.plot(E, dos2_d, label= str(Nth_atom) + '_d')

table = np.column_stack((E, dos1_s, dos1_p, dos1_d, dos2_s, dos2_p, dos2_d))
np.savetxt('LDOS.txt', table, '%.6f', '\t')

plt.legend()
#plt.axis([-7, 5, -3, 10])
plt.xlabel('Energy (eV)')
plt.ylabel('LDOS (State / atom / eV)')
plt.savefig('LDOS.png')
#plt.show()
