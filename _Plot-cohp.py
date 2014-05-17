#!/usr/bin/env python
import sys
import numpy as np
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

if len(list[4]) == 10:
    spin_calc = True
elif len(list[4]) == 5:
    spin_calc = False
else:
    raise Exception("Can't read COHPCAR.lobster properly!")

if spin_calc == False:
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
