#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpltools.style
mpltools.style.use('ggplot')
import re
import argparse


def plot_helper():
    plt.axhline(y=0, c='k')
    plt.axvline(x=0, ls='--', c='k')
    if args.axis_range:
        plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2], args.axis_range[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (States / Unit Cell / eV)')
    plt.legend(loc=0,  fontsize='x-small')
    plt.tight_layout()

parser = argparse.ArgumentParser(description='''Plot the local projected density of states, with
                                                consideration of spin-polarization.''')
parser.add_argument('atom1', type=int, help='first atom to plot')
parser.add_argument('atom2', type=int, help='second atom to plot')
parser.add_argument('-a', '--axis-range', type=eval, help='''the x and y range of axis in the form of
            '[Xmin,Xmax,Ymin,Ymax]'. If ISPIN=2, this option specifies the combined spin''')
parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
parser.add_argument('--LORBIT', type=int, help="manually override LORBIT detection")
parser.add_argument('-i', '--DOSCAR', default='DOSCAR', help="the input DOSCAR file name")
parser.add_argument('-o', '--output-prefix', default='LDOS', help="the output files' prefix")
args = parser.parse_args()
atom1 = args.atom1
atom2 = args.atom2
ISPIN = args.ISPIN
LORBIT = args.LORBIT

with open(args.DOSCAR, 'r') as f:
    DOSCAR = f.readlines()
for i in range(len(DOSCAR)):
    DOSCAR[i] = DOSCAR[i].split()

N_steps = int(DOSCAR[5][2])
Ef = float(DOSCAR[5][3])

if args.ISPIN:
    print "Using user specified ISPIN."
else:
    try:
        with open('OUTCAR', 'r') as f:
            for line in f:
                if 'ISPIN' in line:
                    ISPIN = int(line.split()[2])
    except IOError:
        try:
            with open('INCAR', 'r') as f:
                for line in f:
                    m = re.match(r'\s*ISPIN\s*=\s*(\d)\s*', line)
                    if m:
                        ISPIN = int(m.group(1))
        except IOError:
            raise IOError("Can't determine ISPIN! Either manually specify it, or provide OUTCAR or INCAR")

if args.LORBIT:
    print "Using user specified LORBIT."
else:
    try:
        with open('OUTCAR', 'r') as f:
            for line in f:
                if 'LORBIT' in line:
                    LORBIT = int(line.split()[2])
    except IOError:
        try:
            with open('INCAR', 'r') as f:
                for line in f:
                    m = re.match(r'\s*LORBIT\s*=\s*(\d+)\s*', line)
                    if m:
                        LORBIT = int(m.group(1))
        except IOError:
            raise IOError("Can't determine LORBIT! Either manually specify it, or provide OUTCAR or INCAR")

if ISPIN == 2 and (LORBIT == 11 or LORBIT == 1):
    col_names = ['E', 's_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down',
                 'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z2_up', 'd_z2_down',
                 'd_xz_up', 'd_xz_down', 'd_x2y2_up', 'd_x2y2_down']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    # Spin up for both atoms, above and below
    for i in range(1, 18, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 18, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    if args.axis_range:
        plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
    plt.savefig(args.output_prefix + '-spin-up.png')
    plt.close()

    # Spin down for both atoms, above and below
    for i in range(1, 18, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i + 1], label=col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 18, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i + 1])
    plot_helper()
    if args.axis_range:
        plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
    plt.savefig(args.output_prefix + '-spin-down.png')
    plt.close()

    # Spin up + down for both atoms, above and below
    for i in range(1, 18, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i] + DOS_data_1[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 18, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i] - DOS_data_1[:, i + 1])
    plot_helper()
    plt.savefig(args.output_prefix + '-spin-combined.png')
    plt.close()

elif ISPIN == 2 and (LORBIT == 10 or LORBIT == 0):
    col_names = ['E', 's_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    # Spin up for both atoms, above and below
    for i in range(1, 6, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 6, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    if args.axis_range:
        plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
    plt.savefig(args.output_prefix + '-spin-up.png')
    plt.close()

    # Spin down for both atoms, above and below
    for i in range(1, 6, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i + 1], label=col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 6, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i + 1])
    plot_helper()
    if args.axis_range:
        plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
    plt.savefig(args.output_prefix + '-spin-down.png')
    plt.close()

    # Spin up + down for both atoms, above and below
    for i in range(1, 6, 2):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i] + DOS_data_1[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 6, 2):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i] - DOS_data_1[:, i + 1])
    plot_helper()
    plt.savefig(args.output_prefix + '-spin-combined.png')
    plt.close()

elif ISPIN == 1 and (LORBIT == 11 or LORBIT == 1):
    col_names = ['E', 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2', 'd_xz', 'd_x2y2']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    for i in range(1, 10):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 10):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    plt.savefig(args.output_prefix + '.png')
    plt.close()

elif ISPIN == 1 and (LORBIT == 10 or LORBIT == 0):
    col_names = ['E', 's', 'p', 'd']
    DOS_data_1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
    DOS_data_1[:, 0] -= Ef
    DOS_data_2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
    DOS_data_2[:, 0] -= Ef

    for i in range(1, 4):
        plt.plot(DOS_data_1[:, 0], DOS_data_1[:, i], label=col_names[i])
    ax = plt.gca()
    ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
    for i in range(1, 4):
        plt.plot(DOS_data_2[:, 0], -DOS_data_2[:, i])
    plot_helper()
    plt.savefig(args.output_prefix + '.png')
    plt.close()