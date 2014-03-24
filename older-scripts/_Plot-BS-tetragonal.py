#!/usr/local/python/2.7.1/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

with open('DOSCAR','r') as f:
    for i in range(0,6):
        a = f.readline()
Fermi_E = float(a.split()[3])    # Found in DOSCAR, 6th line, 4th number

f = open('EIGENVAL','r')
file_contents_list_separated_by_line_split = []
for line in f:
    file_contents_list_separated_by_line_split.append(line[0:-1].split())
f.close()

N_steps = int(file_contents_list_separated_by_line_split[5][1])   # How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number
N_bands = int(file_contents_list_separated_by_line_split[5][2])   # How many bands are to be drawn? 6th line, 3rd number

# Get the end point x coordinate of each path on a plot
k_end_Z_G = 0.5
k_end_G_M = np.sqrt(0.5**2 + 0.5**2) + k_end_Z_G
k_end_M_A = 0.5 + k_end_G_M
k_end_A_R = 0.5 + k_end_M_A
k_end_R_X = 0.5 + k_end_A_R
k_end_X_G = 0.5 + k_end_R_X
# Get the array of x coordinates
k_Z_G = np.linspace(0, k_end_Z_G, 10)
k_G_M = np.linspace(k_end_Z_G, k_end_G_M, 10)
k_M_A = np.linspace(k_end_G_M, k_end_M_A, 10)
k_A_R = np.linspace(k_end_M_A, k_end_A_R, 10)
k_R_X = np.linspace(k_end_A_R, k_end_R_X, 10)
k_X_G = np.linspace(k_end_R_X, k_end_X_G, 10)

k = np.hstack((k_Z_G, k_G_M, k_M_A, k_A_R, k_R_X, k_X_G))

# Get energy for each band, each kpoint step
E = np.zeros((N_bands, N_steps))
for n_b in range(0, N_bands):
    for n_s in range(0, N_steps):
        E[n_b, n_s] = float(file_contents_list_separated_by_line_split[8 + n_b + (N_bands + 2) * n_s][1])

E = E - (Fermi_E)
# Relative Fermi energy, choosing the valence band top at Gamma point.
#Fermi_E = E[3][10]

# Plot the bands.
for i in range(0, N_bands):
    plt.plot(k, E[i], label='E')

# Set the axis range.
xlim0 = 0
xlim1 = k_end_X_G
ylim0 = -3
ylim1 = 5
plt.axis([0, xlim1, ylim0, ylim1])

plt.plot([k_end_Z_G, k_end_Z_G], [ylim0, ylim1], 'k--')
plt.plot([k_end_G_M, k_end_G_M], [ylim0, ylim1], 'k--')
plt.plot([k_end_M_A, k_end_M_A], [ylim0, ylim1], 'k--')
plt.plot([k_end_A_R, k_end_A_R], [ylim0, ylim1], 'k--')
plt.plot([k_end_R_X, k_end_R_X], [ylim0, ylim1], 'k--')
plt.plot([k_end_X_G, k_end_X_G], [ylim0, ylim1], 'k--')
plt.plot([xlim0, xlim1],[0, 0], 'k--')

plt.text(-0.02, ylim0*1.2, 'Z')
plt.text(k_end_Z_G-0.02, ylim0*1.2, 'G')
plt.text(k_end_G_M-0.02, ylim0*1.2, 'M')
plt.text(k_end_M_A-0.02, ylim0*1.2, 'A')
plt.text(k_end_A_R-0.02, ylim0*1.2, 'R')
plt.text(k_end_R_X-0.02, ylim0*1.2, 'X')
plt.text(k_end_X_G-0.02, ylim0*1.2, 'G')

frame1 = plt.gca()
frame1.axes.get_xaxis().set_ticks([])
plt.ylabel('Energy (eV)')
plt.savefig('BS.png')

#plt.show()
