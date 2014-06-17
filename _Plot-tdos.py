#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
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

with open('OUTCAR', 'r') as f:
    for line in f:
        if 'ISPIN' in line:
            ISPIN = int(line.split()[2])

if ISPIN == 2:
    # E = np.zeros(N_steps)
    # dos_tot_up = np.zeros(N_steps)
    # dos_tot_down = np.zeros(N_steps)
    # dos_int_up = np.zeros(N_steps)
    # dos_int_down = np.zeros(N_steps)
    # for i in range(0, N_steps):
    #     E[i] = float(DOS_list[i + 6][0]) - Ef
    #     dos_tot_up[i] = float(DOS_list[i + 6][1])
    #     dos_tot_down[i] = float(DOS_list[i + 6][2])
    #     dos_int_up[i] = float(DOS_list[i + 6][3])
    #     dos_int_down[i] = float(DOS_list[i + 6][4])
    col_names = ['E', 'total_up', 'total_down', 'integrated_up', 'integrated_down']
    DOS_data = pd.read_csv('DOSCAR', sep=' *', names=col_names, header=None, skiprows=6, nrows=N_steps)
    DOS_data['E'] -= Ef

    plt.plot(DOS_data['E'], DOS_data['total_up'])
    plt.plot(DOS_data['E'], -DOS_data['total_down'])
    plt.axhline(y=0)
    plt.axis([axis_lim[0], axis_lim[1], -axis_lim[3] / 2., axis_lim[3] / 2.])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS-spin-separated.png')
    plt.close()

    plt.plot(DOS_data['E'], DOS_data['total_up'] + DOS_data['total_down'])
    plt.axis([axis_lim[0], axis_lim[1], 0, axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS-spin-total.png')

    DOS_data.to_csv('TDOS.txt', sep='\t', float_format='%12.6f', index=False)

    slice = DOS_data['integrated_up'].iat[abs(DOS_data['E'] - 0.2).argmin()] - \
            DOS_data['integrated_up'].iat[abs(DOS_data['E'] + 0.2).argmin()] + \
            DOS_data['integrated_down'].iat[abs(DOS_data['E'] - 0.2).argmin()] - \
            DOS_data['integrated_down'].iat[abs(DOS_data['E'] + 0.2).argmin()]

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
    np.savetxt('TDOS@Ef.txt', [slice], '%.6f')
