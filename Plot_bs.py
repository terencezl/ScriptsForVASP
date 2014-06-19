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


# Effective mass calculation funcitons.
def find_band_edges(kp_edge, within, E=E):
    # First identify the k-point where band edges are located: kp_edge
    # Examine valence band edge.
    print 'The possible valence bands are', \
        np.where(np.logical_and(E[:, kp_edge] > -within, E[:, kp_edge] < 0))[0]
    # Examine conduction band edge.
    print 'The possible conduction bands are', \
        np.where(np.logical_and(E[:, kp_edge] < within, E[:, kp_edge] > 0))[0]


def get_effective_mass_reduced(band, kp_start, kp_end, E=E, kp_linearized_array=kp_linearized_array):
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


def main(arguments=['-h']):
    parser = argparse.ArgumentParser(description='''Plot the band structure, with consideration of spin-polarization.''')
    parser.add_argument('-a', '--axis-range', type=eval,
                        help="the x and y range of axis in the form of '[Ymin,Ymax]'.")
    parser.add_argument('-e', '--Ef', type=float, help="manually override Fermi energy detection")
    parser.add_argument('--ISPIN', type=int, help="manually override ISPIN detection")
    parser.add_argument('-i', '--EIGENVAL', default='EIGENVAL', help="the input EIGENVAL file name")
    parser.add_argument('-o', '--output-prefix', default='BS', help="the output files' prefix")
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

    with open(args.EIGENVAL, 'r') as f:
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
                kp_end_letter_list = line.replace(
                    'k-points in units of 2pi/SCALE and weight:', '').strip().split('-')
                break
        for kp in range(N_kps):
            kp_list[kp] = f.next().split()[:3]

    kp_section_start_end_pair_array = np.zeros((N_sections, 2, 3))
    for section in range(N_sections):
        kp_section_start_end_pair_array[section] = [kp_list[N_kps_per_section * section],
                                                    kp_list[N_kps_per_section * (section + 1) - 1]]

    # Gerenate the linearized kp_linearized_array as x-axis.
    section_end_point = 0
    kp_end_point_array = np.zeros(N_sections + 1)
    kp_section_linearized_array = np.zeros((N_sections, N_kps_per_section))
    for section, section_coord_pair in enumerate(kp_section_start_end_pair_array):
        section_end_point_next = section_end_point + np.linalg.norm(
            section_coord_pair[1] - section_coord_pair[0])
        kp_end_point_array[section + 1] = section_end_point_next
        kp_section_linearized_array[section] = np.linspace(section_end_point,
                                                           section_end_point_next, N_kps_per_section)
        section_end_point = section_end_point_next

    kp_linearized_array = kp_section_linearized_array.flatten()

    # Get energy for each band, each kpoint step.
    E = np.zeros((N_bands, N_kps))
    for n_b in range(0, N_bands):
        for n_s in range(0, N_kps):
            E[n_b, n_s] = float(EIGENVAL[8 + n_b + (N_bands + 2) * n_s][1])
    E = E - Ef

    # Plot the bands.
    plt.figure(1)
    ax = plt.subplot(111)
    for band in range(N_bands):
        plt.plot(kp_linearized_array, E[band])

    plt.axis([kp_end_point_array[0], kp_end_point_array[-1], args.axis_range[0], args.axis_range[1]])
    for section_end_point in range(len(kp_end_point_array)):
        plt.axvline(kp_end_point_array[section_end_point], ls='--', c='k', alpha=0.5)
    ax.xaxis.set_ticks(kp_end_point_array)
    ax.xaxis.set_ticklabels(kp_end_letter_list)
    plt.axhline(0, ls='--', c='k', alpha=0.5)
    plt.ylabel('Energy (eV)')
    plt.tight_layout()
    plt.savefig(args.output_prefix + '.png')
    if not plt.isinteractive():
        plt.close()


if __name__ == '__main__':
    main(sys.argv[1:])