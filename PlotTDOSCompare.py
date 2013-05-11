#!/usr/bin/python
#import sys
import matplotlib.pyplot as plt
import numpy as np

N_steps = 301

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

# the total DOS and the integrated DOS
Fermi_E1 = float(list1[5][3])
for i in range( 6, 6 + N_steps ):
    E1.append(float(list1[i][0]) - Fermi_E1)
    dos1_tot.append(float(list1[i][1]))
plt.plot(E1, dos1_tot, label="1")

Fermi_E2 = float(list2[5][3])
for i in range( 6, 6 + N_steps ):
    E2.append(float(list2[i][0]) - Fermi_E2)
    dos2_tot.append(float(list2[i][1]))
plt.plot(E2, dos2_tot, label="2")

#D_dos = np.array(dos1_tot) - np.array(dos2_tot)
#plt.plot (E1, D_dos, label='diff')
 
plt.legend()
#plt.xlim([-15,30])
plt.title("Total DOS")
plt.xlabel("Energy (eV)")
plt.ylabel("DOS")
#plt.savefig('dos.png')
plt.show()
