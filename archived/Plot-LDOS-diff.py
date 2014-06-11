#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

N_steps = 301
Nth_atom = int(sys.argv[1])
#Name_atom = 'Si'

f1 = open('DOSCAR1','rU')
list1 = []; E1 = []; dos1_tot = [];
for line in f1:
    list1.append(line[0:-1].split())
f1.close()

f2 = open('DOSCAR2','rU')
list2 = []; E2 = []; dos2_tot = []
for line in f2:
    list2.append(line[0:-1].split())
f2.close()

# The projected DOS
Ef1 = float(list1[5][3])
E1 = []; dos1_s = []; dos1_p = []; dos1_d = []
for n_s in range(0, N_steps):
    E1.append(float(list1[6+(N_steps+1)*Nth_atom+n_s][0]) - Ef1)
    dos1_s.append(float(list1[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos1_p.append(float(list1[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos1_d.append(float(list1[6+(N_steps+1)*Nth_atom+n_s][3]))
#plt.plot(E1, dos1_s, label='Terence_s')
#plt.plot(E1, dos1_p, label='Terence_p')
#plt.plot(E1, dos1_d, label='Terence_d')

Ef2 = float(list2[5][3])
E2 = []; dos2_s = []; dos2_p = []; dos2_d = []
for n_s in range(0, N_steps):
    E2.append(float(list2[6+(N_steps+1)*Nth_atom+n_s][0]) - Ef2)
    dos2_s.append(float(list2[6+(N_steps+1)*Nth_atom+n_s][1]))
    dos2_p.append(float(list2[6+(N_steps+1)*Nth_atom+n_s][2]))
    dos2_d.append(float(list2[6+(N_steps+1)*Nth_atom+n_s][3]))
#plt.plot(E2, dos2_s, label='Jason_s')
#plt.plot(E2, dos2_p, label='Jason_p')
#plt.plot(E2, dos2_d, label='Jason_d')

D_dos_s = np.array(dos1_s) - np.array(dos2_s)
D_dos_p = np.array(dos1_p) - np.array(dos2_p)
D_dos_d = np.array(dos1_d) - np.array(dos2_d)
plt.plot (E1, D_dos_s, label='s_diff')
plt.plot (E1, D_dos_p, label='s_diff')
plt.plot (E1, D_dos_d, label='s_diff')

plt.legend()
plt.title('Local DOS')
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
#plt.savefig('DOS.png')
plt.show()
