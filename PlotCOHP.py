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
import warnings


def plot_helper_figure_assert(args, ISPIN):
    if ISPIN == 2:
        assert args.figure is None or (isinstance(args.figure, list) and len(args.figure) == 2), \
            'The number of figures should be 2!'
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
    plt.ylabel('-pCOHP (Arbituary Unit / Unit Cell / eV)')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.legend(loc=0,  fontsize='small')
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
            The column names of data from COHPCAR.lobster in a list.
        COHP_data: 2D numpy array
            The data from COHPCAR.lobster.
    Usages
    ------
    main(), main('-h') : as if ``PlotCOHP.py -h`` is executed from command line.
    main(string) : as if ``PlotCOHP.py content_of_string`` is executed from command line.
    """
    arguments = arguments.split()
    parser = argparse.ArgumentParser(description='''Plot the -COHP, with consideration of spin-polarization.''')
    parser.add_argument('bond_to_plot', type=int, help='No. of bond to plot')
    parser.add_argument('-a', '--axis-range', type=eval, help='''the x and y range of axis in the form of
                '[Xmin,Xmax,Ymin,Ymax]'. If ISPIN=2, this option specifies the combined spin.''')
    parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
    parser.add_argument('-i', '--input', metavar='COHPCAR.lobster', default='COHPCAR.lobster', help="the input COHPCAR.lobster file name")
    parser.add_argument('-o', '--output-prefix', default='COHP', help="the output files' prefix")
    parser.add_argument('-f', '--figure', type=eval, help='''the figure number one wishes to plot on,
                                    in the form of '[1,2,...]'. Useful in interactive mode.''')
    args = parser.parse_args(arguments)
    n_bond_to_plot = args.bond_to_plot
    ISPIN = args.ISPIN

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

    plot_helper_figure_assert(args, ISPIN)

    with open(args.input, 'r') as f:
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
        plot_helper_figure(args, ISPIN)
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_up_to_plot], label=col_names[col_up_to_plot])
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_down_to_plot], label=col_names[col_down_to_plot])
        plot_helper_settings(args)
        if args.axis_range:
            plt.axis([args.axis_range[0], args.axis_range[1], args.axis_range[2]/2., args.axis_range[3]/2.])
        plt.savefig(args.output_prefix + '-spin-separated.pdf')
        plot_helper_close()

        # Plot the combined COHP
        plot_helper_figure(args, ISPIN)
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_up_to_plot] - COHP_data[:, col_down_to_plot])
        plot_helper_settings(args)
        plt.savefig(args.output_prefix + '-spin-combined.pdf')
        plot_helper_close()

    elif ISPIN == 1:
        col_names = ['E', 'avg', 'avg_integrated']
        for n_bond in range(1, N_bonds + 1):
            col_names.extend(['No.{0}'.format(n_bond), 'No.{0}_integrated'.format(n_bond)])

        COHP_data = np.array(COHPCAR[data_start_line:data_start_line + N_steps], dtype=float)

        col_to_plot = n_bond_to_plot * 2 + 1
        plot_helper_figure(args, ISPIN)
        plt.plot(COHP_data[:, 0], -COHP_data[:, col_to_plot], label=col_names[col_to_plot])
        plot_helper_settings(args)
        plt.savefig(args.output_prefix + '.pdf')
        plot_helper_close()

    return (col_names, COHP_data)


if __name__ == '__main__':
    main(' '.join(i.replace(' ', '') for i in sys.argv[1:]))