#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import pdb

if len(sys.argv) == 2:
    axis_lim = eval(sys.argv[1])
#else:
#    axis_lim = [-20, 10, -15, 15]

with open('COHPCAR.lobster','r') as f:
    list = f.readlines()
for i in range(0, len(list)):
    list[i] = list[i].split()

N_steps = int(list[1][2])

print list[4]
if len(list[5]) == 13:
    spin_calc = True
elif len(list[5]) == 5:
    spin_calc = False
else:
    raise Exception("Can't read COHPCAR.lobster properly!")

if spin_calc == True:
    E = np.zeros(N_steps)
    cohp1 = np.zeros(N_steps)
    cohp2 = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[i+5][0])
        cohp1[i] = float(list[i+5][3]) + float(list[i+5][9])
        cohp2[i] = float(list[i+5][5]) + float(list[i+5][11])

    plt.plot(E, -cohp1, label="M-M")
    plt.plot(E, -cohp2, label="M-N")
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, ls='--', color='k')
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2], axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('-COHP (Arbituary Unit / Unit Cell / eV)')
    plt.savefig('COHP.png')
    plt.close()

elif spin_calc == False:
    E = np.zeros(N_steps)
    cohp = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[i+4][0])
        cohp[i] = float(list[i+4][1])
    plt.plot(E, -cohp)
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, ls='--', color='k')
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2], axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('-COHP (Arbituary Unit / Unit Cell / eV)')
    plt.savefig('COHP.png')
    plt.close()

