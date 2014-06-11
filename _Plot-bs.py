# !/usr/bin/env python
import sys
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import pandas as pd


if len(sys.argv) == 2 and re.match(r'[.*]', sys.argv[1]):
    [ylim0, ylim1] = eval(sys.argv[1])
else:
    ylim0 = -5
    ylim1 = 5

with open('DOSCAR', 'r') as f:
    for i in range(6):
        a = f.readline()
# Fermi energy. Found in DOSCAR, 6th line, 4th number.
Ef = float(a.split()[3])

with open('EIGENVAL', 'r') as f:
    EIGENVAL = f.readlines()
    for i in range(len(EIGENVAL)):
        EIGENVAL[i] = EIGENVAL[i].split()

# How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number.
N_kps = int(EIGENVAL[5][1])
# How many bands are to be drawn? 6th line, 3rd number.
N_bands = int(EIGENVAL[5][2])

with open('KPOINTS', 'r') as f:
    KPOINTS = f.readlines()
# end_letter_list = KPOINTS[0].strip().split('-')
N_kps_per_section = int(KPOINTS[1])
N_sections = N_kps / N_kps_per_section


# Get the start and end point coordinate of each section. From OUTCAR.
kp_list = [''] * N_kps
with open('OUTCAR', 'r') as f:
    for line in f:
        if re.match(r".*k-points in units of 2pi/SCALE and weight:.*", line):
            head_kp_list_in_cart = line.replace(
                "k-points in units of 2pi/SCALE and weight:", '').strip()
            break
    for kp in range(N_kps):
        kp_list[kp] = f.next().split()[:3]

end_letter_list = head_kp_list_in_cart.split('-')

kp_section_start_end_pair_array = np.zeros((N_sections, 2, 3))
for section in range(6):
    kp_section_start_end_pair_array[section] = [kp_list[N_kps_per_section * section],
                                                kp_list[N_kps_per_section * (section + 1) - 1]]

# Intermediate approach. No longer used.
# kp_section_start_end_array = np.zeros((N_sections * 2, 3))
# for section in range(6):
# kp_section_start_end_array[2 * section] = kp_list[N_kps_per_section * section]
#     kp_section_start_end_array[2 * section + 1] = \
#                                      kp_list[N_kps_per_section * (section + 1) - 1]


# From KPOINTS. No longer used.
# kp_section_start_end_array = []
# for line in KPOINTS[4:]:
#     if line != '\n':
#         kp_section_start_end_array.append(line.split())
#
# kp_section_start_end_pair_array = []
# for i in range(0, len(kp_section_start_end_array) - 2, 2):
#     kp_section_start_end_pair_array.append((np.array(kp_section_start_end_array[i],
#              dtype=float), np.array(kp_section_start_end_array[i+1], dtype=float)))


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
# Relative Fermi energy, choosing the valence band top at Gamma point.
#Ef = E[3][10]


# Plot the bands.
ax = plt.subplot(111)
for band in range(0, N_bands):
    plt.plot(kp_linearized_array, E[band])

# Set the axis properties.
plt.axis([kp_end_point_array[0], kp_end_point_array[-1], ylim0, ylim1])

for section_end_point in range(len(kp_end_point_array)):
    plt.axvline(kp_end_point_array[section_end_point], ls='--', c='k')
    plt.text(kp_end_point_array[section_end_point], -abs(ylim0) * 1.1,
             end_letter_list[section_end_point], ha='center')
plt.axhline(0, ls='--', c='k')

ax.get_xaxis().set_ticks([])
plt.ylabel('Energy (eV)')
plt.savefig('BS.png')




# Effective mass calculation funcitons.
def find_band_edges(kp_edge, within):
    # First identify the k-point where band edges are located: kp_edge
    # Examine valence band edge.
    print "The possible valence bands are", \
        np.where(np.logical_and(E[:, kp_edge] > -within, E[:, kp_edge] < 0))[0]
    # Examine conduction band edge.
    print "The possible conduction bands are", \
        np.where(np.logical_and(E[:, kp_edge] < within, E[:, kp_edge] > 0))[0]


def effective_mass_reduced(band, kp_start, kp_end):
    h_bar = 1.054571726e-34
    e = 1.6021176462e-19
    m_e = 9.10938291e-31
    scaling_const = 6.3743775177e-10

    df = pd.DataFrame({'E': E[band, :]}, index=pd.Series(kp_linearized_array, name='kp'))

    # Decide on the fitting range, characterized by indices.
    p = np.poly1d(np.polyfit(df.index[kp_start:kp_end + 1], df.E.iloc[kp_start:kp_end + 1], 2))
    k_fit = np.linspace(df.index[kp_start], df.index[kp_end], 200)
    plt.plot(k_fit, p(k_fit), lw=2)

    d2E_dk2 = e * p[2] / (2 * np.pi / scaling_const) ** 2
    effective_mass_reduced = (h_bar) ** 2 / d2E_dk2 / m_e
    return effective_mass_reduced
