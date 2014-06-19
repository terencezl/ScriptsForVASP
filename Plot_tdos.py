#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpltools.style
mpltools.style.use('ggplot')
import re
import argparse
import warnings


def plot_helper(args):
    plt.axvline(x=0, ls='--', c='k')
    if args.axis_range:
        plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2], args.axis_range[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.legend(loc=0,  fontsize='small')
    plt.tight_layout()


def main(arguments=['-h']):
    parser = argparse.ArgumentParser(description='''Plot the total density of states, with
                consideration of spin-polarization.''')
    parser.add_argument('-a', '--axis-range', type=eval, help='''the x and y range of axis in the form of
                '[Xmin,Xmax,Ymin,Ymax]'. If ISPIN=2, this option specifies the combined spin.''')
    parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
    parser.add_argument('-i', '--DOSCAR', default='DOSCAR', help="the input DOSCAR file name")
    parser.add_argument('-o', '--output-prefix', default='TDOS', help="the output files' prefix")
    args = parser.parse_args(arguments)
    ISPIN = args.ISPIN

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

    if ISPIN == 2:
        col_names = ['E', 'total_up', 'total_down', 'integrated_up', 'integrated_down']
        DOS_data = np.array(DOSCAR[6:6+N_steps], dtype=float)
        DOS_data[:, 0] -= Ef

        # Plot the separated TDOS
        plt.figure(1)
        plt.plot(DOS_data[:, 0], DOS_data[:, 1], label='spin up')
        plt.plot(DOS_data[:, 0], DOS_data[:, 2], label='spin down')
        plot_helper(args)
        if args.axis_range:
            plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2], args.axis_range[3] / 2.])
        plt.savefig(args.output_prefix + '-spin-separated.png')
        if not plt.isinteractive():
            plt.close()

        # Plot the combined TDOS
        plt.figure(2)
        plt.plot(DOS_data[:, 0], DOS_data[:, 1] + DOS_data[:, 2], label='spin up + down')
        plot_helper(args)
        plt.savefig(args.output_prefix + '-spin-combined.png')
        if not plt.isinteractive():
            plt.close()

        np.savetxt(args.output_prefix + '.txt', DOS_data, '%15.6E', header=' '.join(col_names))
        energy_slice = DOS_data[abs(DOS_data[0] - 0.2).argmin(), 3] - \
                DOS_data[abs(DOS_data[0] + 0.2).argmin(), 3] + \
                DOS_data[abs(DOS_data[0] - 0.2).argmin(), 4] - \
                DOS_data[abs(DOS_data[0] + 0.2).argmin(), 4]
        np.savetxt(args.output_prefix + '@Ef.txt', [energy_slice], '%15.6E')

    elif ISPIN == 1:
        col_names = ['E', 'total', 'integrated']
        DOS_data = np.array(DOSCAR[6:6+N_steps], dtype=float)
        DOS_data[:, 0] -= Ef

        plt.figure(1)
        plt.plot(DOS_data[:, 0], DOS_data[:, 1])
        plot_helper(args)
        plt.savefig(args.output_prefix + '.png')
        if not plt.isinteractive():
            plt.close()

        np.savetxt(args.output_prefix + '.txt', DOS_data, '%15.6E', header=' '.join(col_names))
        energy_slice = DOS_data[abs(DOS_data[0] - 0.2).argmin(), 2] - \
                       DOS_data[abs(DOS_data[0] + 0.2).argmin(), 2]
        np.savetxt(args.output_prefix + '@Ef.txt', [energy_slice], '%15.6E')


if __name__ == '__main__':
    main(sys.argv[1:])