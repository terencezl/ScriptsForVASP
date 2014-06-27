#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
if not matplotlib.is_interactive():
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpltools.style
mpltools.style.use('ggplot')
import re
import argparse


# Effective mass calculation funcitons.
def find_band_edges(kp_edge, prec_range, E):
    """
    Given the k-point index number on the x-axis, search for the indices of bands
    right below or above Ef.
    Best used in an interactive Python interpreter and having pyplot.ion().

    Parameters
    ----------
    kp_edge : int
        The k-point index number that corresponds to the band edges (VBM and CBM).
    prec_range : float
        The searching range in eV.
    E : 2D numpy array
        The 2D array that contains eigenvalues. Each row denotes a band, and each column a k-point.
    """
    # Examine valence band edge.
    print 'The possible valence bands are', \
        np.where(np.logical_and(E[:, kp_edge] > -prec_range, E[:, kp_edge] < 0))[0]
    # Examine conduction band edge.
    print 'The possible conduction bands are', \
        np.where(np.logical_and(E[:, kp_edge] < prec_range, E[:, kp_edge] > 0))[0]


def get_effective_mass_reduced(band, kp_start, kp_end, kp_linearized_array, E):
    """
    Given the band index number, k-point start and end indices, fit the included curve
    to a 2nd-order polynomial, and obtain the effective mass of the carrier electron or hole.
    Best used in an interactive Python interpreter and having pyplot.ion().

    Parameters
    ----------
    band : int
        The band index number of interest.
    kp_start : int
        The index number of the starting k-point.
    kp_end : int
        The index number of the ending k-point.
    kp_linearized_array : 1D numpy array
        The full x-axis of k-points.
    E : 2D numpy array
        The 2D array that contains eigenvalues. Each row denotes a band, and each column a k-point.

    Returns
    -------
    effective_mass_reduced : float
        The reduced effective mass.
    """
    h_bar = 1.054571726e-34
    e = 1.6021176462e-19
    m_e = 9.10938291e-31
    scaling_const = 6.3743775177e-10

    # Decide on the fitting range, characterized by indices.
    selected_kp_array = kp_linearized_array[kp_start:kp_end + 1]
    selected_energy_array = E[band, kp_start:kp_end + 1]
    p = np.poly1d(np.polyfit(selected_kp_array, selected_energy_array, 2))
    axis_fitted = -p[1]/2/p[2]
    axis_actual = selected_kp_array[selected_energy_array.argmin() if p[2] > 0 else selected_energy_array.argmax()]
    print "The fitted x coord at energy extrema is {0}, and the actual is {1}.".format(axis_fitted, axis_actual)
    k_fit = np.linspace(kp_linearized_array[kp_start], kp_linearized_array[kp_end], 200)
    plt.plot(k_fit, p(k_fit), lw=2)

    d2E_dk2 = e * p[2] / (2 * np.pi / scaling_const) ** 2
    effective_mass_reduced = h_bar ** 2 / d2E_dk2 / m_e
    return effective_mass_reduced


def plot_helper_figure(args, ISPIN):
    if ISPIN == 2:
        assert args.figure is None or (isinstance(args.figure, list) and len(args.figure) == 2), \
            'The number of figures should be 2!'
    elif ISPIN == 1:
        assert args.figure is None or (isinstance(args.figure, list) and len(args.figure) == 1), \
            'The number of figures should be 1!'
    if args.figure is None:
        plt.figure()
    else:
        plt.figure(args.figure.pop(0))


def plot_helper_close():
    if not matplotlib.is_interactive():
        plt.close()


def main(arguments='-h'):
    """
    The main function that does the data grepping and plotting.

    :param arguments: string
        Command line arguments and options. Typically ' '.join(sys.argv[1:]).
    :return:
        kp_end_symbol_list: list
            The symbols of points of interest in the 1st BZ, as a sequence of x-axis labels.
        kp_end_point_list: list
            The positional values of the ending k-point of each line-mode section on x-axis.
        kp_linearized_array: 1D numpy array
            The full x-axis of k-points.
        E: 2D numpy array
            The 2D array that contains eigenvalues. Each row denotes a band, and each column a k-point.
         
    Usages
    ------
    main(), main('-h') : as if ``PlotBS.py -h`` is executed from command line.
    main(string) : as if ``PlotBS.py content_of_string`` is executed from command line.
    """
    arguments = arguments.split()
    parser = argparse.ArgumentParser(description='''Plot the band structure, with consideration of spin-polarization.''')
    parser.add_argument('-a', '--axis-range', type=eval,
                        help="the x and y range of axis in the form of '[Ymin,Ymax]'.")
    parser.add_argument('-e', '--Ef', type=float, help="manually override Fermi energy detection")
    parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
    parser.add_argument('-i', '--input', metavar='EIGENVAL', default='EIGENVAL', help="the input EIGENVAL file name")
    parser.add_argument('-o', '--output-prefix', default='BS', help="the output files' prefix")
    parser.add_argument('-f', '--figure', type=eval, help='''the figure number one wishes to plot on,
                                    in the form of '[1,2,...]'. Useful in interactive mode.''')
    args = parser.parse_args(arguments)
    Ef = args.Ef
    ISPIN = args.ISPIN

    if args.Ef:
        print "Using user specified Ef."
    else:
        try:
            with open('DOSCAR', 'r') as f:
                for i in range(6):
                    line = f.readline()
            # Fermi energy. Found in DOSCAR, 6th line, 4th number.
            Ef = float(line.split()[3])
        except IOError:
            raise IOError("Can't determine Ef! Either manually specify it, or provide DOSCAR")

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

    with open(args.input, 'r') as f:
        EIGENVAL = f.readlines()
    for i in range(len(EIGENVAL)):
        EIGENVAL[i] = EIGENVAL[i].split()

    # How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number.
    N_kps = int(EIGENVAL[5][1])
    # How many bands are to be drawn? 6th line, 3rd number.
    N_bands = int(EIGENVAL[5][2])

    with open('KPOINTS', 'r') as f:
        KPOINTS = f.readlines()
    N_kps_per_section = int(KPOINTS[1])
    N_sections = N_kps / N_kps_per_section

    # Get the start and end point coordinate of each section. From OUTCAR.
    kp_list = np.zeros((N_kps, 3))
    with open('OUTCAR', 'r') as f:
        for line in f:
            if re.match(r".*k-points in units of 2pi/SCALE and weight:.*", line):
                kp_end_symbol_list = line.replace(
                    'k-points in units of 2pi/SCALE and weight:', '').strip().split('-')
                break
        for kp in range(N_kps):
            kp_list[kp] = f.next().split()[:3]

    kp_section_start_end_pair_array = np.zeros((N_sections, 2, 3))
    for section in range(N_sections):
        kp_section_start_end_pair_array[section] = [kp_list[N_kps_per_section * section],
                                                    kp_list[N_kps_per_section * (section + 1) - 1]]

    # Gerenate the linearized kp_linearized_array as x-axis.
    kp_end_point = 0
    kp_end_point_list = [0] * (N_sections + 1)
    kp_section_linearized_array = np.zeros((N_sections, N_kps_per_section))
    for section, section_coord_pair in enumerate(kp_section_start_end_pair_array):
        kp_end_point_next = kp_end_point + np.linalg.norm(
            section_coord_pair[1] - section_coord_pair[0])
        kp_end_point_list[section + 1] = kp_end_point_next
        kp_section_linearized_array[section] = np.linspace(kp_end_point,
                                                           kp_end_point_next, N_kps_per_section)
        kp_end_point = kp_end_point_next

    kp_linearized_array = kp_section_linearized_array.flatten()


    # Get energy for each band, each kpoint step.
    E = np.zeros((N_bands, N_kps))
    for n_b in range(0, N_bands):
        for n_s in range(0, N_kps):
            E[n_b, n_s] = float(EIGENVAL[8 + n_b + (N_bands + 2) * n_s][1])
    E -= Ef

    # Plot the bands.
    plot_helper_figure(args, ISPIN)
    ax = plt.subplot(111)
    for band in range(N_bands):
        plt.plot(kp_linearized_array, E[band])

    plt.xlim(kp_end_point_list[0], kp_end_point_list[-1])
    ax.xaxis.set_ticks(kp_end_point_list)
    ax.xaxis.set_ticklabels(kp_end_symbol_list)
    plt.axhline(0, ls='--', c='k', alpha=0.5)
    if args.axis_range:
        plt.ylim(args.axis_range[0], args.axis_range[1])
    for kp_end_point in range(len(kp_end_point_list)):
        plt.axvline(kp_end_point_list[kp_end_point], ls='--', c='k', alpha=0.5)
    plt.ylabel('Energy (eV)')
    try:
        plt.tight_layout()
    except RuntimeError:
        print "Tight layout failed... Not a big deal though."
    plt.savefig(args.output_prefix + '.png')
    plot_helper_close()

    return kp_end_symbol_list, kp_end_point_list, kp_linearized_array, E


if __name__ == '__main__':
    main(' '.join(sys.argv[1:]))