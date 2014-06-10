#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import pdb

if len(sys.argv) == 2:
    axis_lim = eval(sys.argv[1])
else:
    axis_lim = [-25, 10, 0, 15]

with open('DOSCAR','r') as f:
    list = f.readlines()
for i in range(0, len(list)):
    list[i] = list[i].split()

N_steps = int(list[5][2])
Fermi_E = float(list[5][3])

if len(list[6]) == 3:
    is_spin = False
elif len(list[6]) == 5:
    is_spin = True
else:
    raise Exception("Can't read DOSCAR properly!")

if is_spin:
    E = np.zeros(N_steps)
    dos_tot_up = np.zeros(N_steps)
    dos_tot_down = np.zeros(N_steps)
    dos_int_up = np.zeros(N_steps)
    dos_int_down = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[i+6][0]) - Fermi_E
        dos_tot_up[i] = float(list[i+6][1])
        dos_tot_down[i] = float(list[i+6][2])
        dos_int_up[i] = float(list[i+6][3])
        dos_int_down[i] = float(list[i+6][4])
    plt.plot(E, dos_tot_up)
    plt.plot(E, -dos_tot_down)
    plt.axhline(y=0)
    plt.axis([axis_lim[0], axis_lim[1], -axis_lim[3]/2., axis_lim[3]/2.])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS-spin-separated.png')
    plt.close()

    plt.plot(E, dos_tot_up + dos_tot_down)
    plt.axis([axis_lim[0], axis_lim[1], 0, axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS-spin-total.png')

    table = np.column_stack((E, dos_tot_up + dos_tot_down, np.zeros(len(E))))
    np.savetxt('TDOS.txt', table, '%.6f', '\t')
    slice = dos_int_up[np.abs(np.array(E) - 0.2).argmin()] - dos_int_up[np.abs(np.array(E) + 0.2).argmin()] + dos_int_down[np.abs(np.array(E) - 0.2).argmin()] - dos_int_down[np.abs(np.array(E) + 0.2).argmin()]
    #slice = dos_tot_up[np.abs(np.array(E)).argmin()] + dos_tot_down[np.abs(np.array(E)).argmin()]
    np.savetxt('TDOS@Ef.txt', [slice], '%.6f')

else:
    E = np.zeros(N_steps)
    dos_tot = np.zeros(N_steps)
    dos_int = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[i+6][0]) - Fermi_E
        dos_tot[i] = float(list[i+6][1])
        dos_int[i] = float(list[i+6][2])
    plt.plot(E, dos_tot)
    plt.axis([axis_lim[0], axis_lim[1], 0, axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS.png')

    table = np.column_stack((E, dos_tot))
    np.savetxt('TDOS.txt', table, '%.6f', '\t')
    slice = dos_int[np.abs(np.array(E) - 0.2).argmin()] - dos_int[np.abs(np.array(E) + 0.2).argmin()]
    #slice = dos_tot[np.abs(np.array(E)).argmin()]
    np.savetxt('TDOS@Ef.txt', [slice], '%.6f')
