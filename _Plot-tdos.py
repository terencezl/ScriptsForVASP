#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re


def plot_helper():
    plt.axvline(x=0, ls='--', c='k')
    plt.axis([axis_lim[0], axis_lim[1], 0, axis_lim[2]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.legend(loc=0,  fontsize='x-small')
    plt.tight_layout()


if len(sys.argv) == 2 and re.match(r'\[.*\]', sys.argv[1]):
    axis_lim = eval(sys.argv[1])
else:
    axis_lim = [-25, 10, 10]

with open('DOSCAR', 'r') as f:
    DOSCAR = f.readlines()
for i in range(len(DOSCAR)):
    DOSCAR[i] = DOSCAR[i].split()

N_steps = int(DOSCAR[5][2])
Ef = float(DOSCAR[5][3])

with open('OUTCAR', 'r') as f:
    for line in f:
        if 'ISPIN' in line:
            ISPIN = int(line.split()[2])

if ISPIN == 2:
    col_names = ['E', 'total_up', 'total_down', 'integrated_up', 'integrated_down']
    DOS_data = np.array(DOSCAR[6:6+N_steps], dtype=float)
    DOS_data[:, 0] -= Ef

    # Plot the separated TDOS
    plt.plot(DOS_data[:, 0], DOS_data[:, 1], label='spin up')
    plt.plot(DOS_data[:, 0], DOS_data[:, 2], label='spin down')
    plot_helper()
    plt.savefig('TDOS-spin-separated.png')
    plt.close()

    # Plot the combined TDOS
    plt.plot(DOS_data[:, 0], DOS_data[:, 1] + DOS_data[:, 2])
    plot_helper()
    plt.savefig('TDOS-spin-combined.png')
    plt.close()

    np.savetxt('TDOS.txt', DOS_data, '%15.6E', header=' '.join(col_names))
    energy_slice = DOS_data[abs(DOS_data[0] - 0.2).argmin(), 3] - \
            DOS_data[abs(DOS_data[0] + 0.2).argmin(), 3] + \
            DOS_data[abs(DOS_data[0] - 0.2).argmin(), 4] - \
            DOS_data[abs(DOS_data[0] + 0.2).argmin(), 4]
    np.savetxt('TDOS@Ef.txt', [energy_slice], '%15.6E')

else:
    col_names = ['E', 'total', 'integrated']
    DOS_data = np.array(DOSCAR[6:6+N_steps], dtype=float)
    DOS_data[:, 0] -= Ef

    plt.plot(DOS_data[:, 0], DOS_data[:, 1])
    plot_helper()
    plt.savefig('TDOS.png')
    plt.close()

    np.savetxt('TDOS.txt', DOS_data, '%15.6E', header=' '.join(col_names))
    energy_slice = DOS_data[abs(DOS_data[0] - 0.2).argmin(), 2] - \
                   DOS_data[abs(DOS_data[0] + 0.2).argmin(), 2]
    np.savetxt('TDOS@Ef.txt', [energy_slice], '%15.6E')