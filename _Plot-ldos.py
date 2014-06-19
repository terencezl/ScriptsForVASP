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
    plt.axis([axis_lim[0], axis_lim[1], -axis_lim[2], axis_lim[2]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (States / Unit Cell / eV)')
    plt.legend(loc=0, fontsize='x-small')
    plt.tight_layout()


if re.match(r'\[.*\]', sys.argv[1]) and re.match(r'\[.*\]', sys.argv[2]):
    axis_lim = eval(sys.argv[1])
    atom_1, atom_2 = eval(sys.argv[2])
else:
    axis_lim = [-20, 10, 2]
    atom_1, atom_2 = 1, 2

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
        if 'LORBIT' in line:
            LORBIT = int(line.split()[2])

if ISPIN == 2 and (LORBIT == 11 or LORBIT == 1):
    col_names = ['E', 's_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down',
                 'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z2_up', 'd_z2_down',
                 'd_xz_up', 'd_xz_down', 'd_x2y2_up', 'd_x2y2_down']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_1):(6 + (N_steps + 1) * atom_1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_2):(6 + (N_steps + 1) * atom_2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    # Spin up for both atoms, above and below
    for i in range(1, 18, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 18, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    plt.savefig('LDOS-spin-up.png')
    plt.close()

    # Spin down for both atoms, above and below
    for i in range(1, 18, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i + 1], label=col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 18, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i + 1])
    plot_helper()
    plt.savefig('LDOS-spin-down.png')
    plt.close()

    # Spin up + down for both atoms, above and below
    for i in range(1, 18, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i] + DOS_data_1[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 18, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i] - DOS_data_1[:, i + 1])
    plot_helper()
    plt.savefig('LDOS-spin-combined.png')
    plt.close()

elif ISPIN == 2 and (LORBIT == 10 or LORBIT == 0):
    col_names = ['E', 's_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_1):(6 + (N_steps + 1) * atom_1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_2):(6 + (N_steps + 1) * atom_2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    # Spin up for both atoms, above and below
    for i in range(1, 6, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 6, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    plt.savefig('LDOS-spin-up.png')
    plt.close()

    # Spin down for both atoms, above and below
    for i in range(1, 6, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i + 1], label=col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 6, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i + 1])
    plot_helper()
    plt.savefig('LDOS-spin-down.png')
    plt.close()

    # Spin up + down for both atoms, above and below
    for i in range(1, 6, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i] + DOS_data_1[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 6, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i] - DOS_data_1[:, i + 1])
    plot_helper()
    plt.savefig('LDOS-spin-combined.png')
    plt.close()

elif ISPIN == 1 and (LORBIT == 11 or LORBIT == 1):
    col_names = ['E', 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2', 'd_xz', 'd_x2y2']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_1):(6 + (N_steps + 1) * atom_1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_2):(6 + (N_steps + 1) * atom_2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    for i in range(1, 10):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 10):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    plt.savefig('LDOS.png')
    plt.close()

elif ISPIN == 1 and (LORBIT == 10 or LORBIT == 0):
    col_names = ['E', 's', 'p', 'd']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_1):(6 + (N_steps + 1) * atom_1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom_2):(6 + (N_steps + 1) * atom_2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    for i in range(1, 4):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 4):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    plt.savefig('LDOS.png')
    plt.close()