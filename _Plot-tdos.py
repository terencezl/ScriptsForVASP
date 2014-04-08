#!/usr/local/python/2.7.1/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
#import pdb

if len(sys.argv) == 2:
    axis_lim = eval(sys.argv[1])
else:
    axis_lim = [-20, 10, 0, 15]

with open('DOSCAR','r') as f:
    list = f.readlines()
for i in range(0, len(list)):
    list[i] = list[i].split()

N_steps = int(list[5][2])
Fermi_E = float(list[5][3])

if len(list[6]) == 3:
    spin_calc = False
elif len(list[6]) == 5:
    spin_calc = True
else:
    print "Can't read DOSCAR properly!"
    sys.exit(0)

# The total DOS and the integrated DOS

if spin_calc == True:
    E = np.zeros(N_steps)
    dos_tot_up = np.zeros(N_steps)
    dos_tot_down = np.zeros(N_steps)
    dos_int_up = np.zeros(N_steps)
    dos_int_down = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[i+6][0]) - Fermi_E
        dos_tot_up[i] = float(list[i+6][1])
        dos_tot_down[i] = float(list[i+6][2])
        dos_int_up[i] = float(list[i+6][3])
        dos_int_down[i] = float(list[i+6][4])
    plt.plot(E, dos_tot_up)
    plt.plot(E, -dos_tot_down)
    plt.axhline(y=0)
    plt.axis([axis_lim[0], axis_lim[1], -axis_lim[3]/2., axis_lim[3]/2.])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS-spin-separated.png')
    plt.close()
    
    plt.plot(E, dos_tot_up + dos_tot_down)
    plt.axis([axis_lim[0], axis_lim[1], 0, axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS-spin-total.png')

    table = np.column_stack((E, dos_tot_up, dos_int_down, dos_tot_up + dos_tot_down))
    np.savetxt('TDOS.txt', table, '%.6f', '\t')
    #slice = dos_int[np.abs(np.array(E) - 0.2).argmin()] - dos_int[np.abs(np.array(E) + 0.2).argmin()]
    #slice = dos_tot[np.abs(np.array(E)).argmin()]
    #np.savetxt('TDOS@Ef.txt', [slice], '%.6f')
    
else:
    E = np.zeros(N_steps)
    dos_tot = np.zeros(N_steps)
    dos_int = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[i+6][0]) - Fermi_E
        dos_tot[i] = float(list[i+6][1])
        dos_int[i] = float(list[i+6][2])
    plt.plot(E, dos_tot)
    #plt.plot(E, dos_int, label="Integrated")
    plt.axis([axis_lim[0], axis_lim[1], 0, axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('TDOS (States / Unit Cell / eV)')
    plt.savefig('TDOS.png')

    table = np.column_stack((E, dos_tot))
    np.savetxt('TDOS.txt', table, '%.6f', '\t')
    #slice = dos_int[np.abs(np.array(E) - 0.2).argmin()] - dos_int[np.abs(np.array(E) + 0.2).argmin()]
    #slice = dos_tot[np.abs(np.array(E)).argmin()]
    #np.savetxt('TDOS@Ef.txt', [slice], '%.6f')
