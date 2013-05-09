#!/usr/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

# How many KPs in total? Can be found in EIGENVAL, 6th line, 2nd number
N_steps = 10
# How many bands are to be drawn? 6th line, 3rd number
N_bands = 8
# Found in DOSCAR, 6th line, 4th number
# Fermi_E1 = 5.49111961 
# Control group
# Fermi_E2 = 5.35618101

f1 = open('EIGENVAL1','rU')
list1 = []
for line in f1:
    list1.append(line[0:-1].split())
f1.close()

# Control group
f2 = open('EIGENVAL2','rU')
list2 = []
for line in f2:
    list2.append(line[0:-1].split())
f2.close()

#KPOINTS
#Si L-G-X-K-G
#Cartesian
#0.5     0.5     0.5
#0.0     0.0     0.0

#0.0     0.0     0.0
#0.0     1.0     0.0

#0.0     1.0     0.0
#0.75    0.75    0.0

#0.75    0.75    0.0
#0.0     0.0     0.0

# Get the end point x coordinate of each path on a plot
k_end_L_G = np.sqrt(0.5**2 + 0.5**2 + 0.5**2)
k_end_G_X = np.sqrt(1) + k_end_L_G
k_end_X_K = np.sqrt(0.75**2 + 0.25**2) + k_end_G_X
k_end_K_G = np.sqrt(0.75**2 + 0.75**2) + k_end_X_K
# Get the array of x coordinates
k_L_G = np.linspace(0, k_end_L_G, 10)
k_G_X = np.linspace(k_end_L_G, k_end_G_X, 10)
k_X_K = np.linspace(k_end_G_X, k_end_X_K, 10)
k_K_G = np.linspace(k_end_X_K, k_end_K_G, 10)

k1 = np.hstack((k_L_G, k_G_X, k_X_K, k_K_G))

# Control group
k_end_M_K = np.sqrt(0.25**2 + 0.25**2) + k_end_G_X
k_M_K = np.linspace(k_end_G_X, k_end_M_K, 10)

k2 = np.hstack((k_L_G, k_G_X, k_M_K, k_K_G))


E1 = []
for i in range(0, N_bands):
    E1.append([])

E2 = []
for i in range(0, N_bands):
    E2.append([])

# Get every E points for N_bands bands at N * N_steps KPs
for n_s in range(0, N_steps * 4):
    for n_b in range(0, N_bands):
        E1[n_b].append(float(list1[8+n_b+(N_bands+2)*n_s][1]))
Fermi_E1 = E1[3][10]
E1 = np.array(E1) - Fermi_E1

for n_s in range(0, N_steps * 4):
    for n_b in range(0, N_bands):
        E2[n_b].append(float(list2[8+n_b+(N_bands+2)*n_s][1]))
Fermi_E2 = E2[3][10]
E2 = np.array(E2) - Fermi_E2

dE = E1 - E2

# Plot the bands
for i in range(0, N_bands):
    h1, = plt.plot(k1, E1[i], 'k-', label='E1_'+str(i+1))

for i in range(0, N_bands):
    h2, = plt.plot(k2, E2[i], 'r*', label='E2_'+str(i+1))

for i in range(0, N_bands):
    plt.plot(k2[:20], dE[i][:20], 'b-', label='dE'+str(i+1))
    h3, = plt.plot(k2[30:], dE[i][30:], 'b-', label='dE'+str(i+1))

plt.legend([h1, h2, h3], ["Terence", "Jason", "Comparison"])
# Set the axis range
xlim0 = 0
xlim1 = k_end_K_G
ylim0 = -15
ylim1 = 15
plt.axis([0, xlim1, ylim0, ylim1])

plt.plot([k_end_L_G, k_end_L_G], [ylim0, ylim1], 'k--')
plt.plot([k_end_G_X, k_end_G_X], [ylim0, ylim1], 'k--')
plt.plot([k_end_X_K, k_end_X_K], [ylim0, ylim1], 'k--')
plt.plot([k_end_K_G, k_end_K_G], [ylim0, ylim1], 'k--')
plt.plot([xlim0, xlim1],[0, 0], 'k--')

plt.text(-0.02, ylim0-1, 'L')
plt.text(k_end_L_G-0.02, ylim0-1, 'G')
plt.text(k_end_G_X-0.02, ylim0-1, 'X')
plt.text(k_end_X_K-0.02, ylim0-1, 'K')
plt.text(k_end_K_G-0.02, ylim0-1, 'G')
plt.text(xlim1-0.3, 0.3, r'$E_F$')

plt.text(k_end_G_X-0.02, ylim0-2, 'M')


frame1 = plt.gca()
frame1.axes.get_xaxis().set_ticks([])
plt.title("Band Structure")
plt.ylabel('Energy (eV)')
plt.savefig('bandstructure.png')
plt.show()
