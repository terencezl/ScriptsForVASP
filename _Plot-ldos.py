#!/usr/local/python/2.7.1/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) == 3:
    atom_1st = int(sys.argv[1])
    atom_2nd = int(sys.argv[2])
else:
    atom_1st = 1
    atom_2nd = 2

if len(sys.argv) == 2:
    axis_lim = eval(sys.argv[1])
else:
    axis_lim = [-20, 10, -10, 10]

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

if spin_calc == True:
    E = np.zeros(N_steps)
    dos1_up_s = np.zeros(N_steps)
    dos1_up_p = np.zeros(N_steps)
    dos1_up_d = np.zeros(N_steps)
    dos1_down_s = np.zeros(N_steps)
    dos1_down_p = np.zeros(N_steps)
    dos1_down_d = np.zeros(N_steps)
    dos2_up_s = np.zeros(N_steps)
    dos2_up_p = np.zeros(N_steps)
    dos2_up_d = np.zeros(N_steps)
    dos2_down_s = np.zeros(N_steps)
    dos2_down_p = np.zeros(N_steps)
    dos2_down_d = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[6+(N_steps+1)*atom_1st+i][0]) - Fermi_E
        dos1_up_s[i] = float(list[6+(N_steps+1)*atom_1st+i][1])
        dos1_up_p[i] = float(list[6+(N_steps+1)*atom_1st+i][3])
        dos1_up_d[i] = float(list[6+(N_steps+1)*atom_1st+i][5])
        dos1_down_s[i] = float(list[6+(N_steps+1)*atom_1st+i][2])
        dos1_down_p[i] = float(list[6+(N_steps+1)*atom_1st+i][4])
        dos1_down_d[i] = float(list[6+(N_steps+1)*atom_1st+i][6])
        dos2_up_s[i] = float(list[6+(N_steps+1)*atom_2nd+i][1])
        dos2_up_p[i] = float(list[6+(N_steps+1)*atom_2nd+i][3])
        dos2_up_d[i] = float(list[6+(N_steps+1)*atom_2nd+i][5])
        dos2_down_s[i] = float(list[6+(N_steps+1)*atom_2nd+i][2])
        dos2_down_p[i] = float(list[6+(N_steps+1)*atom_2nd+i][4])
        dos2_down_d[i] = float(list[6+(N_steps+1)*atom_2nd+i][6])

    # spin-up
    plt.plot(E, dos1_up_s, label= str(atom_1st) + '_s')
    plt.plot(E, dos1_up_p, label= str(atom_1st) + '_p')
    plt.plot(E, dos1_up_d, label= str(atom_1st) + '_d')
    plt.plot(E, -dos2_up_s, label= str(atom_2nd) + '_s')
    plt.plot(E, -dos2_up_p, label= str(atom_2nd) + '_p')
    plt.plot(E, -dos2_up_d, label= str(atom_2nd) + '_d')
    plt.legend(loc=0)
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2]/2., axis_lim[3]/2.])
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (State / atom / eV)')
    plt.savefig('LDOS-spin-up.png')
    plt.close()

    # spin-down
    plt.plot(E, dos1_down_s, label= str(atom_1st) + '_s')
    plt.plot(E, dos1_down_p, label= str(atom_1st) + '_p')
    plt.plot(E, dos1_down_d, label= str(atom_1st) + '_d')
    plt.plot(E, -dos2_down_s, label= str(atom_2nd) + '_s')
    plt.plot(E, -dos2_down_p, label= str(atom_2nd) + '_p')
    plt.plot(E, -dos2_down_d, label= str(atom_2nd) + '_d')
    plt.legend(loc=0)
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2]/2., axis_lim[3]/2.])
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (State / atom / eV)')
    plt.savefig('LDOS-spin-down.png')
    plt.close()

    # spin-total
    plt.plot(E, dos1_up_s + dos1_down_s, label= str(atom_1st) + '_s')
    plt.plot(E, dos1_up_p + dos1_down_p, label= str(atom_1st) + '_p')
    plt.plot(E, dos1_up_d + dos1_down_d, label= str(atom_1st) + '_d')
    plt.plot(E, -dos2_up_s - dos2_down_s, label= str(atom_2nd) + '_s')
    plt.plot(E, -dos2_up_p - dos2_down_p, label= str(atom_2nd) + '_p')
    plt.plot(E, -dos2_up_d - dos2_down_d, label= str(atom_2nd) + '_d')
    plt.legend(loc=0)
    plt.axis(axis_lim)
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (State / atom / eV)')
    plt.savefig('LDOS-spin-total.png')
    plt.close()

#    table = np.column_stack((E, dos1_s, dos1_p, dos1_d, dos2_s, dos2_p, dos2_d))
#    np.savetxt('LDOS.txt', table, '%.6f', '\t')

else:
    E = np.zeros(N_steps)
    dos1_s = np.zeros(N_steps)
    dos1_p = np.zeros(N_steps)
    dos1_d = np.zeros(N_steps)
    dos2_s = np.zeros(N_steps)
    dos2_p = np.zeros(N_steps)
    dos2_d = np.zeros(N_steps)
    for i in range(0, N_steps):
        E[i] = float(list[6+(N_steps+1)*atom_1st+i][0]) - Fermi_E
        dos1_s[i] = float(list[6+(N_steps+1)*atom_1st+i][1])
        dos1_p[i] = float(list[6+(N_steps+1)*atom_1st+i][2])
        dos1_d[i] = float(list[6+(N_steps+1)*atom_1st+i][3])
        dos2_s[i] = float(list[6+(N_steps+1)*atom_2nd+i][1])
        dos2_p[i] = float(list[6+(N_steps+1)*atom_2nd+i][2])
        dos2_d[i] = float(list[6+(N_steps+1)*atom_2nd+i][3])

    plt.plot(E, dos1_s, label= str(atom_1st) + '_s')
    plt.plot(E, dos1_p, label= str(atom_1st) + '_p')
    plt.plot(E, dos1_d, label= str(atom_1st) + '_d')
    plt.plot(E, -dos2_s, label= str(atom_2nd) + '_s')
    plt.plot(E, -dos2_p, label= str(atom_2nd) + '_p')
    plt.plot(E, -dos2_d, label= str(atom_2nd) + '_d')
    plt.legend(loc=0)
    plt.axis(axis_lim)
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (State / atom / eV)')
    plt.savefig('LDOS-spin-up.png')
    plt.close()

#    table = np.column_stack((E, dos1_s, dos1_p, dos1_d, dos2_s, dos2_p, dos2_d))
#    np.savetxt('LDOS.txt', table, '%.6f', '\t')
