#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpltools.style
mpltools.style.use('ggplot')
import re


def plot_helper():
    plt.axhline(y=0, c='k')
    plt.axvline(x=0, ls='--', c='k')
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2], axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('-pCOHP (Arbituary Unit / Unit Cell / eV)')
    plt.legend(loc=0,  fontsize='small')
    plt.tight_layout()


if re.match(r'\[.*\]', sys.argv[1]) and re.match(r'\d', sys.argv[2]):
    axis_lim = eval(sys.argv[1])
    n_bond_to_plot = int(sys.argv[2])
    filename = sys.argv[3]
    ISPIN = int(sys.argv[4])
else:
    axis_lim = [-25, 10, -2, 5]
    n_bond_to_plot = 1
    filename = 'COHPCAR.lobster'

with open(filename, 'r') as f:
    COHPCAR = f.readlines()

for N_headerlines, line in enumerate(COHPCAR):
    if re.match(r'No\.\d*:.*\(.*\)', line):
        break
for N_bonds, line in enumerate(COHPCAR[N_headerlines:]):
    if not re.match(r'No\.\d*:.*\(.*\)', line):
        break
data_start_line = N_headerlines + N_bonds

for i in range(len(COHPCAR)):
    COHPCAR[i] = COHPCAR[i].split()

N_steps = int(COHPCAR[1][2])

if 'ISPIN' not in globals():
    try:
        with open('OUTCAR', 'r') as f:
            for line in f:
                if 'ISPIN' in line:
                    ISPIN = int(line.split()[2])
    except IOError:
        try:
            with open('INCAR', 'r') as f:
                for line in f:
                    if 'ISPIN' in line:
                        ISPIN = int(line.split()[-1])
        except IOError:
            raise IOError('No ISPIN value determined!')

if ISPIN == 2:
    col_names = ['E', 'avg_up', 'avg_integrated_up']
    for n_bond in range(1, N_bonds + 1):
        col_names.extend(['No.{0}_up'.format(n_bond), 'No.{0}_integrated_up'.format(n_bond)])
    col_names.extend(['avg_down', 'avg_integrated_down'])
    for n_bond in range(1, N_bonds + 1):
        col_names.extend(['No.{0}_down'.format(n_bond), 'No.{0}_integrated_down'.format(n_bond)])
    COHP_data = np.array(COHPCAR[data_start_line:data_start_line + N_steps], dtype=float)

    col_up_to_plot = n_bond_to_plot * 2 + 1
    col_down_to_plot = (n_bond_to_plot + N_bonds + 1) * 2 + 1

    # Plot the separated COHP
    plt.plot(COHP_data[:, 0], -COHP_data[:, col_up_to_plot], label=col_names[col_up_to_plot])
    plt.plot(COHP_data[:, 0], -COHP_data[:, col_down_to_plot], label=col_names[col_down_to_plot])
    plot_helper()
    plt.savefig(filename + '-spin-separated.png')
    plt.close()

    # Plot the combined COHP
    plt.plot(COHP_data[:, 0], -COHP_data[:, col_up_to_plot] - COHP_data[:, col_down_to_plot])
    plot_helper()
    plt.savefig(filename + '-spin-combined.png')
    plt.close()

elif ISPIN == 1:
    col_names = ['E', 'avg', 'avg_integrated']
    for n_bond in range(1, N_bonds + 1):
        col_names.extend(['No.{0}'.format(n_bond), 'No.{0}_integrated'.format(n_bond)])

    COHP_data = np.array(COHPCAR[data_start_line:data_start_line + N_steps], dtype=float)

    col_to_plot = n_bond_to_plot * 2 + 1
    plt.plot(COHP_data[:, 0], -COHP_data[:, col_to_plot], label=col_names[col_to_plot])
    plot_helper()
    plt.savefig(filename + '.png')
    plt.close()