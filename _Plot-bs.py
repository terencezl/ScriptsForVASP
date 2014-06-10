#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

if len(sys.argv) == 2:
    [ylim0, ylim1] = eval(sys.argv[1])
else:
    ylim0 = -5
    ylim1 = 5


with open('DOSCAR','r') as f:
    for i in range(6):
        a = f.readline()
# Fermi energy. Found in DOSCAR, 6th line, 4th number.
Ef = float(a.split()[3])


with open('EIGENVAL','r') as f:
    EIGENVAL = f.readlines()
    for i in range(len(EIGENVAL)):
        EIGENVAL[i] = EIGENVAL[i].split()

# How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number.
N_kps = int(EIGENVAL[5][1])
# How many bands are to be drawn? 6th line, 3rd number.
N_bands = int(EIGENVAL[5][2])


with open('KPOINTS','r') as f:
    KPOINTS = f.readlines()
# end_letter_list = KPOINTS[0].strip().split('-')
N_kps_per_section = int(KPOINTS[1])
N_sections = N_kps/N_kps_per_section


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
    kp_section_start_end_pair_array[section] = [ kp_list[N_kps_per_section * section], kp_list[N_kps_per_section * (section + 1) - 1] ]

# Intermediate approach. No longer used.
# kp_section_start_end_array = np.zeros((N_sections * 2, 3))
# for section in range(6):
#     kp_section_start_end_array[2 * section] = kp_list[N_kps_per_section * section]
#     kp_section_start_end_array[2 * section + 1] = kp_list[N_kps_per_section * (section + 1) - 1]


# From KPOINTS. No longer used.
# kp_section_start_end_array = []
# for line in KPOINTS[4:]:
#     if line != '\n':
#         kp_section_start_end_array.append(line.split())
#
# kp_section_start_end_pair_array = []
# for i in range(0, len(kp_section_start_end_array) - 2, 2):
#     kp_section_start_end_pair_array.append((np.array(kp_section_start_end_array[i], dtype=float), np.array(kp_section_start_end_array[i+1], dtype=float)))


# Gerenate the linearized kp_array as x-axis.
section_end_point = 0
kp_end_point_array = np.zeros(N_sections + 1)
kp_section_linearized_array = np.zeros((N_sections, N_kps_per_section))
for section, section_coord_pair in enumerate(kp_section_start_end_pair_array):
    section_end_point_next = section_end_point + np.linalg.norm(section_coord_pair[1] - section_coord_pair[0])
    kp_end_point_array[section + 1] = section_end_point_next
    kp_section_linearized_array[section] = np.linspace(section_end_point, section_end_point_next, N_kps_per_section)
    section_end_point = section_end_point_next

kp_array = kp_section_linearized_array.flatten()


# Get energy for each band, each kpoint step.
E = np.zeros((N_bands, N_kps))
for n_b in range(0, N_bands):
    for n_s in range(0, N_kps):
        E[n_b, n_s] = float(EIGENVAL[8 + n_b + (N_bands + 2) * n_s][1])

E = E - Ef
# Relative Fermi energy, choosing the valence band top at Gamma point.
#Ef = E[3][10]


# Plot the bands.
for i in range(0, N_bands):
    plt.plot(kp_array, E[i])

# Set the axis properties.
plt.axis([kp_end_point_array[0], kp_end_point_array[-1], ylim0, ylim1])

for section_end_point in range(len(kp_end_point_array)):
    plt.axvline(kp_end_point_array[section_end_point], ls='--', c='k')
    plt.text(kp_end_point_array[section_end_point] - 0.02, -abs(ylim0) * 1.1, end_letter_list[section_end_point])
plt.axhline(0, ls='--', c='k')

ax = plt.gca()
ax.get_xaxis().set_ticks([])
plt.ylabel('Energy (eV)')
plt.savefig('BS.png')