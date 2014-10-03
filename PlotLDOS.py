#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
if not matplotlib.is_interactive():
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:
    plt.style.use('ggplot')
except AttributeError:
    print "If you upgrade to matplotlib 1.4 and I will change the style to ggplot, just prettier."
import re
import argparse


def plot_helper_figure_assert(args, ISPIN):
    if ISPIN == 2:
        assert args.figure is None or (isinstance(args.figure, list) and len(args.figure) == 3), \
            'The number of figures should be 3!'
    elif ISPIN == 1:
        assert args.figure is None or (isinstance(args.figure, list) and len(args.figure) == 1), \
            'The number of figures should be 1!'


def plot_helper_figure(args):
    if args.figure is None:
        plt.figure()
    else:
        plt.figure(args.figure.pop(0))
        

def plot_helper_settings(args):
    plt.axhline(y=0, c='k')
    plt.axvline(x=0, ls='--', c='k')
    if args.axis_range:
        plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2], args.axis_range[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (States / Unit Cell / eV)')
    plt.legend(loc=0,  fontsize='x-small')
    try:
        plt.tight_layout()
    except RuntimeError:
        print "Tight layout failed... Not a big deal though."


def plot_helper_close():
    if not matplotlib.is_interactive():
        plt.close()


def main(arguments='-h'):
    """
    The main function that does the data grepping and plotting.

    :param arguments: string
        Command line arguments and options. Typically ' '.join(sys.argv[1:]).
    :return:
        col_names: list
            The column names of data from DOSCAR in a list.
        DOS_data1: 2D numpy array
            The data of the 1st atom from DOSCAR.
        DOS_data2: 2D numpy array
            The data of the 2nd atom from DOSCAR.
    Usages
    ------
    main(), main('-h') : as if ``PlotLDOS.py -h`` is executed from command line.
    main(string) : as if ``PlotLDOS.py content_of_string`` is executed from command line.
    """
    arguments = arguments.split()
    parser = argparse.ArgumentParser(description='''Plot the local projected density of states, with
                                                    consideration of spin-polarization.''')
    parser.add_argument('atom1', metavar='ATOM1', type=int, help='first atom to plot')
    parser.add_argument('atom2', metavar='ATOM2', type=int, help='second atom to plot')
    parser.add_argument('-a', '--axis-range', type=eval, help='''the x and y range of axis in the form of
                '[Xmin,Xmax,Ymin,Ymax]'. If ISPIN=2, this option specifies the combined spin''')
    parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
    parser.add_argument('--LORBIT', type=int, help="manually override LORBIT detection")
    parser.add_argument('-i', '--input', metavar='DOSCAR', default='DOSCAR', help="the input DOSCAR file name")
    parser.add_argument('-o', '--output-prefix', default='LDOS', help="the output files' prefix")
    parser.add_argument('-f', '--figure', type=eval, help='''the figure number one wishes to plot on,
                                    in the form of '[1,2,...]'. Useful in interactive mode.''')
    args = parser.parse_args(arguments)
    atom1 = args.atom1
    atom2 = args.atom2
    ISPIN = args.ISPIN
    LORBIT = args.LORBIT

    with open(args.input, 'r') as f:
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

    plot_helper_figure_assert(args, ISPIN)

    if ISPIN == 2 and (LORBIT == 11 or LORBIT == 1):
        col_names = ['E', 's_up', 's_down', 'p_y_up', 'p_y_down', 'p_z_up', 'p_z_down', 'p_x_up', 'p_x_down',
                     'd_xy_up', 'd_xy_down', 'd_yz_up', 'd_yz_down', 'd_z2_up', 'd_z2_down',
                     'd_xz_up', 'd_xz_down', 'd_x2y2_up', 'd_x2y2_down']
        DOS_data1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
        DOS_data1[:, 0] -= Ef
        DOS_data2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
        DOS_data2[:, 0] -= Ef

        # Spin up for both atoms, above and below
        plot_helper_figure(args)
        for i in range(1, 18, 2):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i], label=col_names[i])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 18, 2):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i])
        plot_helper_settings(args)
        if args.axis_range:
            plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
        plt.savefig(args.output_prefix + '-spin-up.pdf')
        plot_helper_close()

        # Spin down for both atoms, above and below
        plot_helper_figure(args)
        for i in range(1, 18, 2):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i + 1], label=col_names[i + 1])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 18, 2):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i + 1])
        plot_helper_settings(args)
        if args.axis_range:
            plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
        plt.savefig(args.output_prefix + '-spin-down.pdf')
        plot_helper_close()

        # Spin up + down for both atoms, above and below
        plot_helper_figure(args)
        for i in range(1, 18, 2):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i] + DOS_data1[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 18, 2):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i] - DOS_data2[:, i + 1])
        plot_helper_settings(args)
        plt.savefig(args.output_prefix + '-spin-combined.pdf')
        plot_helper_close()

    elif ISPIN == 2 and (LORBIT == 10 or LORBIT == 0):
        col_names = ['E', 's_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']
        DOS_data1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
        DOS_data1[:, 0] -= Ef
        DOS_data2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
        DOS_data2[:, 0] -= Ef

        # Spin up for both atoms, above and below
        plot_helper_figure(args)
        for i in range(1, 6, 2):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i], label=col_names[i])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 6, 2):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i])
        plot_helper_settings(args)
        if args.axis_range:
            plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
        plt.savefig(args.output_prefix + '-spin-up.pdf')
        plot_helper_close()

        # Spin down for both atoms, above and below
        plot_helper_figure(args)
        for i in range(1, 6, 2):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i + 1], label=col_names[i + 1])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 6, 2):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i + 1])
        plot_helper_settings(args)
        if args.axis_range:
            plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
        plt.savefig(args.output_prefix + '-spin-down.pdf')
        plot_helper_close()

        # Spin up + down for both atoms, above and below
        plot_helper_figure(args)
        for i in range(1, 6, 2):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i] + DOS_data1[:, i + 1], label=col_names[i] + '+' + col_names[i + 1])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 6, 2):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i] - DOS_data2[:, i + 1])
        plot_helper_settings(args)
        plt.savefig(args.output_prefix + '-spin-combined.pdf')
        plot_helper_close()

    elif ISPIN == 1 and (LORBIT == 11 or LORBIT == 1):
        col_names = ['E', 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2', 'd_xz', 'd_x2y2']
        DOS_data1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
        DOS_data1[:, 0] -= Ef
        DOS_data2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
        DOS_data2[:, 0] -= Ef

        plot_helper_figure(args)
        for i in range(1, 10):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i], label=col_names[i])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 10):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i])
        plot_helper_settings(args)
        plt.savefig(args.output_prefix + '.pdf')
        plot_helper_close()

    elif ISPIN == 1 and (LORBIT == 10 or LORBIT == 0):
        col_names = ['E', 's', 'p', 'd']
        DOS_data1 = np.array(DOSCAR[(6 + (N_steps + 1) * atom1):(6 + (N_steps + 1) * atom1 + N_steps)], dtype=float)
        DOS_data1[:, 0] -= Ef
        DOS_data2 = np.array(DOSCAR[(6 + (N_steps + 1) * atom2):(6 + (N_steps + 1) * atom2 + N_steps)], dtype=float)
        DOS_data2[:, 0] -= Ef

        plot_helper_figure(args)
        for i in range(1, 4):
            plt.plot(DOS_data1[:, 0], DOS_data1[:, i], label=col_names[i])
        ax = plt.gca()
        ax.set_color_cycle(plt.rcParams['axes.color_cycle'])
        for i in range(1, 4):
            plt.plot(DOS_data2[:, 0], -DOS_data2[:, i])
        plot_helper_settings(args)
        plt.savefig(args.output_prefix + '.pdf')
        plot_helper_close()

    return col_names, DOS_data1, DOS_data2


if __name__ == '__main__':
    main(' '.join(i.replace(' ', '') for i in sys.argv[1:]))