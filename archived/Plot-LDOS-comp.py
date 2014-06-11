#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

N_steps = 301
fname = sys.argv[1]
metal = sys.argv[2]
cryst_struct = sys.argv[3]
f = open(fname,'rU')
list = []
for line in f:
    list.append(line[0:-1].split())
f.close()

# The projected DOS
Nth_atom = 1
Ef = float(list[5][3])
E = []; dos_s = []; dos_p = []; dos_d = []
for n_s in range(0, N_steps):
    E.append(float(list[6+(N_steps+1)*Nth_atom+n_s][0]) - Ef)
    dos_s.append(float(list[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos_p.append(float(list[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos_d.append(float(list[6+(N_steps+1)*Nth_atom+n_s][3]))
#plt.plot(E, dos_s, label= metal+'_s')
plt.plot(E, dos_p, label= metal+'_p')
plt.plot(E, dos_d, label= metal+'_d')

Nth_atom = 2
Ef = float(list[5][3])
E = []; dos_s = []; dos_p = []; dos_d = []
for n_s in range(0, N_steps):
    E.append(float(list[6+(N_steps+1)*Nth_atom+n_s][0]) - Ef)
    dos_s.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos_p.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos_d.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][3]))
plt.plot(E, dos_s, label='N_s')
plt.plot(E, dos_p, label='N_p')
#plt.plot(E, dos_d, label='N_d')

plt.legend()
#plt.xlim([-17,7])
plt.ylim([-3, 5])
#plt.title('Local DOS')
plt.xlabel('Energy (eV)')
plt.ylabel('LDOS (State/atom/eV)')
plt.savefig('LDOS-'+metal+'N-'+cryst_struct+'.png')
#plt.show()
