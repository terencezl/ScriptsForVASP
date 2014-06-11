#!/usr/bin/env python
import sys
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import pdb

if len(sys.argv) == 2:
    axis_lim = eval(sys.argv[1])
# else:
#    axis_lim = [-20, 10, -15, 15]

with open('COHPCAR.lobster', 'r') as f:
    COHP_list = f.readlines()
for i in range(0, len(COHP_list)):
    COHP_list[i] = COHP_list[i].split()

N_steps = int(COHP_list[1][2])

print COHP_list[4]
if len(COHP_list[5]) == 13:
    is_spin = True
elif len(COHP_list[5]) == 5:
    is_spin = False
else:
    raise Exception("Can't read COHPCAR.lobster properly!")

if is_spin:
    E = np.zeros(N_steps)
    cohp1 = np.zeros(N_steps)
    cohp2 = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(COHP_list[i + 5][0])
        cohp1[i] = float(COHP_list[i + 5][3]) + float(COHP_list[i + 5][9])
        cohp2[i] = float(COHP_list[i + 5][5]) + float(COHP_list[i + 5][11])

    plt.plot(E, -cohp1, label="M-M")
    plt.plot(E, -cohp2, label="M-N")
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, ls='--', color='k')
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2], axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('-COHP (Arbituary Unit / Unit Cell / eV)')
    plt.savefig('COHP.png')
    plt.close()

elif not is_spin:
    E = np.zeros(N_steps)
    cohp = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(COHP_list[i + 4][0])
        cohp[i] = float(COHP_list[i + 4][1])
    plt.plot(E, -cohp)
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, ls='--', color='k')
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2], axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('-COHP (Arbituary Unit / Unit Cell / eV)')
    plt.savefig('COHP.png')
    plt.close()

