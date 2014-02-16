#!/usr/local/python/2.7.1/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

N_steps = 1501
fname = sys.argv[1]
metal = sys.argv[3]
cryst_struct = sys.argv[2]
strain_type = sys.argv[4]
strain_ratio = sys.argv[5]

f = open(fname,'rU')
list = []
for line in f:
    list.append(line[0:-1].split())
f.close()

# The projected DOS
Nth_atom = 1
Fermi_E = float(list[5][3])
E = []; dos1_s = []; dos1_p = []; dos1_d = []
for n_s in range(0, N_steps):
    E.append(float(list[6+(N_steps+1)*Nth_atom+n_s][0]) - Fermi_E)
    dos1_s.append(float(list[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos1_p.append(float(list[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos1_d.append(float(list[6+(N_steps+1)*Nth_atom+n_s][3]))
plt.plot(E, dos1_s, label= metal+'_s')
plt.plot(E, dos1_p, label= metal+'_p')
plt.plot(E, dos1_d, label= metal+'_d')

Nth_atom = 5
#E = [];
dos2_s = []; dos2_p = []; dos2_d = []
for n_s in range(0, N_steps):
#    E2.append(float(list[6+(N_steps+1)*Nth_atom+n_s][0]) - Fermi_E)
    dos2_s.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos2_p.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos2_d.append(-float(list[6+(N_steps+1)*Nth_atom+n_s][3]))
plt.plot(E, dos2_s, label='N_s')
plt.plot(E, dos2_p, label='N_p')
plt.plot(E, dos2_d, label='N_d')

table = np.column_stack((E, dos1_s, dos1_p, dos1_d, dos2_s, dos2_p, dos2_d))
np.savetxt('LDOS-'+metal+'N-'+cryst_struct+'-'+strain_type+'-'+strain_ratio+'.txt', table, '%.6f', '\t')

#plt.legend()
plt.axis([-7, 5, -3, 10])
#plt.title('Local DOS')
plt.xlabel('Energy (eV)')
plt.ylabel('LDOS (State / atom / eV)')
plt.savefig('LDOS-'+metal+'N-'+cryst_struct+'-'+strain_type+'-'+strain_ratio+'.png')
#plt.show()
