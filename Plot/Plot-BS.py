#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

N_steps = 10   # How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number
N_bands = 8    # How many bands are to be drawn? 6th line, 3rd number
Fermi_E = 5.49111961  # Found in DOSCAR, 6th line, 4th number

f = open('EIGENVAL','rU')
list = []
for line in f:
    list.append(line[0:-1].split())
f.close()

# Get the end point x coordinate of each path on a plot
k_end_L_G = np.sqrt(0.5**2 + 0.5**2 + 0.5**2)
k_end_G_X = np.sqrt(1) + k_end_L_G
k_end_M_K = np.sqrt(0.25**2 + 0.25**2) + k_end_G_X
k_end_K_G = np.sqrt(0.75**2 + 0.75**2) + k_end_M_K
# Get the array of x coordinates
k_L_G = np.linspace(0, k_end_L_G, 10)
k_G_X = np.linspace(k_end_L_G, k_end_G_X, 10)
k_M_K = np.linspace(k_end_G_X, k_end_M_K, 10)
k_K_G = np.linspace(k_end_M_K, k_end_K_G, 10)

k = np.hstack((k_L_G, k_G_X, k_M_K, k_K_G))

E = []
for i in range(0, N_bands):
    E.append([])

# Get every E points for N_bands bands at N * N_steps KPs.
for n_s in range(0, N_steps * 4):
    for n_b in range(0, N_bands):
        E[n_b].append(float(list[8+n_b+(N_bands+2)*n_s][1]))
# Relative Fermi energy, choosing the valence band top at Gamma point.
Fermi_E = E[3][10]
E = np.array(E) - Fermi_E

# Plot the bands.
for i in range(0, N_bands):
    plt.plot(k, E[i], label='E')

# Set the axis range.
xlim0 = 0
xlim1 = k_end_K_G
ylim0 = -15
ylim1 = 12
plt.axis([0, xlim1, ylim0, ylim1])

plt.plot([k_end_L_G, k_end_L_G], [ylim0, ylim1], 'k--')
plt.plot([k_end_G_X, k_end_G_X], [ylim0, ylim1], 'k--')
plt.plot([k_end_M_K, k_end_M_K], [ylim0, ylim1], 'k--')
plt.plot([k_end_K_G, k_end_K_G], [ylim0, ylim1], 'k--')
plt.plot([xlim0, xlim1],[0, 0], 'k--')

plt.text(-0.02, ylim0-1, 'L')
plt.text(k_end_L_G-0.02, ylim0-1, 'G')
plt.text(k_end_G_X-0.02, ylim0-1, 'X')
plt.text(k_end_M_K-0.02, ylim0-1, 'U,K')
plt.text(k_end_K_G-0.02, ylim0-1, 'G')
#plt.text(xlim1-0.3, 0.3, r'$E_F$')

frame1 = plt.gca()
frame1.axes.get_xaxis().set_ticks([])
plt.title("Si Band Structure")
plt.ylabel('Energy (eV)')
#plt.savefig('GaAs_BS.png')
plt.show()
