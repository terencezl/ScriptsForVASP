#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) == 2:
    [ylim0, ylim1] = eval(sys.argv[1])
else:
    ylim0 = -5
    ylim1 = 5

with open('DOSCAR','r') as f:
    for i in range(0,6):
        a = f.readline()
Fermi_E = float(a.split()[3])    # Found in DOSCAR, 6th line, 4th number

with open('EIGENVAL','r') as f:
    EIGENVAL = f.readlines()
    for i in range(len(EIGENVAL)):
        EIGENVAL[i] = EIGENVAL[i].split()

N_steps = int(EIGENVAL[5][1])   # How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number
N_bands = int(EIGENVAL[5][2])   # How many bands are to be drawn? 6th line, 3rd number

with open('KPOINTS','r') as f:
    KPOINTS = f.readlines()

end_letter_list = KPOINTS[0][:-1].split('-')

# get the end point x coordinate of each path on a plot
k_start_and_end = []
for line in KPOINTS[4:]:
    if line != '\n':
        k_start_and_end.append(line.split())

k_start_end_tuple_list = []
for i in range(0, len(k_start_and_end) - 1, 2):
    k_start_end_tuple_list.append((np.array(k_start_and_end[i], dtype=float), np.array(k_start_and_end[i+1], dtype=float)))

point_per_path = int(N_steps) / len(k_start_end_tuple_list)

end_point = 0
end_point_list = [end_point]
k_list_sections = []
for i in k_start_end_tuple_list:
    end_point_next = end_point + np.linalg.norm(i[1] - i[0])
    end_point_list.append(end_point_next)
    k_list_sections.append(np.linspace(end_point, end_point_next, point_per_path))
    end_point = end_point_next

k_list = np.array(k_list_sections).flatten()

# get energy for each band, each kpoint step
E = np.zeros((N_bands, N_steps))
for n_b in range(0, N_bands):
    for n_s in range(0, N_steps):
        E[n_b, n_s] = float(EIGENVAL[8 + n_b + (N_bands + 2) * n_s][1])

E = E - (Fermi_E)
# Relative Fermi energy, choosing the valence band top at Gamma point.
#Fermi_E = E[3][10]

# plot the bands
for i in range(0, N_bands):
    plt.plot(k_list, E[i])

# set the axis properties
plt.axis([end_point_list[0], end_point_list[-1], ylim0, ylim1])

for i in range(len(end_point_list)):
    plt.plot([end_point_list[i], end_point_list[i]], [ylim0, ylim1], 'k--')
    plt.text(end_point_list[i] - 0.02, ylim0 * 1.1, end_letter_list[i])
plt.plot([end_point_list[0], end_point_list[-1]], [0, 0], 'k--')

frame = plt.gca()
frame.axes.get_xaxis().set_ticks([])
plt.ylabel('Energy (eV)')
plt.savefig('BS.png')

#plt.show()
