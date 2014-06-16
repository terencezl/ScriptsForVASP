#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

if len(sys.argv) == 2 and re.match(r'\[.*\]', sys.argv[1]):
    axis_lim = eval(sys.argv[1])
else:
    axis_lim = [-25, 10, 0, 15]

with open('DOSCAR', 'r') as f:
    DOS_list = f.readlines()
for i in range(len(DOS_list)):
    DOS_list[i] = DOS_list[i].split()

N_steps = int(DOS_list[5][2])
Ef = float(DOS_list[5][3])

if len(DOS_list[6]) == 3:
    is_spin = False
elif len(DOS_list[6]) == 5:
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
        E[i] = float(DOS_list[i + 6][0]) - Ef
        dos_tot_up[i] = float(DOS_list[i + 6][1])
        dos_tot_down[i] = float(DOS_list[i + 6][2])
        dos_int_up[i] = float(DOS_list[i + 6][3])
        dos_int_down[i] = float(DOS_list[i + 6][4])
    plt.plot(E, dos_tot_up)
    plt.plot(E, -dos_tot_down)
    plt.axhline(y=0)
    plt.axis([axis_lim[0], axis_lim[1], -axis_lim[3] / 2., axis_lim[3] / 2.])
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
    slice = dos_int_up[np.abs(np.array(E) - 0.2).argmin()] - dos_int_up[np.abs(np.array(E) + 0.2).argmin()] + \
            dos_int_down[np.abs(np.array(E) - 0.2).argmin()] - dos_int_down[np.abs(np.array(E) + 0.2).argmin()]
    # slice = dos_tot_up[np.abs(np.array(E)).argmin()] + dos_tot_down[np.abs(np.array(E)).argmin()]
    np.savetxt('TDOS@Ef.txt', [slice], '%.6f')

else:
    E = np.zeros(N_steps)
    dos_tot = np.zeros(N_steps)
    dos_int = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(DOS_list[i + 6][0]) - Ef
        dos_tot[i] = float(DOS_list[i + 6][1])
        dos_int[i] = float(DOS_list[i + 6][2])
    plt.plot(E, dos_tot)
    plt.axis([axis_lim[0], axis_lim[1], 0, axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS.png')

    table = np.column_stack((E, dos_tot))
    np.savetxt('TDOS.txt', table, '%.6f', '\t')
    slice = dos_int[np.abs(np.array(E) - 0.2).argmin()] - dos_int[np.abs(np.array(E) + 0.2).argmin()]
    # slice = dos_tot[np.abs(np.array(E)).argmin()]
    np.savetxt('TDOS@Ef.txt', [slice], '%.6f')
