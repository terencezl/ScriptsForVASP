#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) == 3:
    atom_1st = int(sys.argv[1])
    atom_2nd = int(sys.argv[2])
else:
    atom_1st = 1
    atom_2nd = 5

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

if len(list[6]) == 5:
    is_spin = True
    if len(list[6+N_steps+1]) == 19:
        is_m_decomposed = True
    elif len(list[6+N_steps+1]) == 7:
        is_m_decomposed = False
    else:
        raise Exception("Can't read DOSCAR properly!")
elif len(list[6]) == 3:
    is_spin = False
    if len(list[6+N_steps+1]) == 10:
        is_m_decomposed = True
    elif len(list[6+N_steps+1]) == 4:
        is_m_decomposed = False
    else:
        raise Exception("Can't read DOSCAR properly!")
else:
    raise Exception("Can't read DOSCAR properly!")

if is_spin and is_m_decomposed:
    E = np.zeros(N_steps)
    dos1_s = np.zeros(N_steps)
    dos1_p_y = np.zeros(N_steps)
    dos1_p_z = np.zeros(N_steps)
    dos1_p_x = np.zeros(N_steps)
    dos1_d_xy = np.zeros(N_steps)
    dos1_d_yz = np.zeros(N_steps)
    dos1_d_z2 = np.zeros(N_steps)
    dos1_d_xz = np.zeros(N_steps)
    dos1_d_x2_y2 = np.zeros(N_steps)
    dos2_s = np.zeros(N_steps)
    dos2_p_y = np.zeros(N_steps)
    dos2_p_z = np.zeros(N_steps)
    dos2_p_x = np.zeros(N_steps)
    dos2_d_xy = np.zeros(N_steps)
    dos2_d_yz = np.zeros(N_steps)
    dos2_d_z2 = np.zeros(N_steps)
    dos2_d_xz = np.zeros(N_steps)
    dos2_d_x2_y2 = np.zeros(N_steps)

    for i in range(0, N_steps):
        E[i] = float(list[6+(N_steps+1)*atom_1st+i][0]) - Fermi_E
        dos1_s[i] = float(list[6+(N_steps+1)*atom_1st+i][1]) + float(list[6+(N_steps+1)*atom_1st+i][2])
        dos1_p_y[i] = float(list[6+(N_steps+1)*atom_1st+i][3]) + float(list[6+(N_steps+1)*atom_1st+i][4])
        dos1_p_z[i] = float(list[6+(N_steps+1)*atom_1st+i][5]) + float(list[6+(N_steps+1)*atom_1st+i][6])
        dos1_p_x[i] = float(list[6+(N_steps+1)*atom_1st+i][7]) + float(list[6+(N_steps+1)*atom_1st+i][8])
        dos1_d_xy[i] = float(list[6+(N_steps+1)*atom_1st+i][9]) + float(list[6+(N_steps+1)*atom_1st+i][10])
        dos1_d_yz[i] = float(list[6+(N_steps+1)*atom_1st+i][11]) + float(list[6+(N_steps+1)*atom_1st+i][12])
        dos1_d_z2[i] = float(list[6+(N_steps+1)*atom_1st+i][13]) + float(list[6+(N_steps+1)*atom_1st+i][14])
        dos1_d_xz[i] = float(list[6+(N_steps+1)*atom_1st+i][15]) + float(list[6+(N_steps+1)*atom_1st+i][16])
        dos1_d_x2_y2[i] = float(list[6+(N_steps+1)*atom_1st+i][17]) + float(list[6+(N_steps+1)*atom_1st+i][18])

        dos2_s[i] = float(list[6+(N_steps+1)*atom_2nd+i][1]) + float(list[6+(N_steps+1)*atom_2nd+i][2])
        dos2_p_y[i] = float(list[6+(N_steps+1)*atom_2nd+i][3]) + float(list[6+(N_steps+1)*atom_2nd+i][4])
        dos2_p_z[i] = float(list[6+(N_steps+1)*atom_2nd+i][5]) + float(list[6+(N_steps+1)*atom_2nd+i][6])
        dos2_p_x[i] = float(list[6+(N_steps+1)*atom_2nd+i][7]) + float(list[6+(N_steps+1)*atom_2nd+i][8])
        dos2_d_xy[i] = float(list[6+(N_steps+1)*atom_2nd+i][9]) + float(list[6+(N_steps+1)*atom_2nd+i][10])
        dos2_d_yz[i] = float(list[6+(N_steps+1)*atom_2nd+i][11]) + float(list[6+(N_steps+1)*atom_2nd+i][12])
        dos2_d_z2[i] = float(list[6+(N_steps+1)*atom_2nd+i][13]) + float(list[6+(N_steps+1)*atom_2nd+i][14])
        dos2_d_xz[i] = float(list[6+(N_steps+1)*atom_2nd+i][15]) + float(list[6+(N_steps+1)*atom_2nd+i][16])
        dos2_d_x2_y2[i] = float(list[6+(N_steps+1)*atom_2nd+i][17]) + float(list[6+(N_steps+1)*atom_2nd+i][18])


    # E = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 0], dtype=float ) - Fermi_E
    # dos1_s = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 1], dtype=float ) + \
    #          np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 2], dtype=float )
    # dos1_p_y = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 3], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 4], dtype=float )
    # dos1_p_z = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 5], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 6], dtype=float )
    # dos1_p_x = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 7], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 8], dtype=float )
    # dos1_d_xy = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 9], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 10], dtype=float )
    # dos1_d_yz = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 11], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 12], dtype=float )
    # dos1_d_z2 = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 13], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 14], dtype=float )
    # dos1_d_xz = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 15], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 16], dtype=float )
    # dos1_d_x2_y2 = np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 17], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_1st):(6+(N_steps+1)*atom_1st+N_steps), 18], dtype=float )

    # dos2_s = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 1], dtype=float ) + \
    #          np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 2], dtype=float )
    # dos2_p_y = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 3], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 4], dtype=float )
    # dos2_p_z = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 5], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 6], dtype=float )
    # dos2_p_x = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 7], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 8], dtype=float )
    # dos2_d_xy = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 9], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 10], dtype=float )
    # dos2_d_yz = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 11], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 12], dtype=float )
    # dos2_d_z2 = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 13], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 14], dtype=float )
    # dos2_d_xz = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 15], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 16], dtype=float )
    # dos2_d_x2_y2 = np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 17], dtype=float ) + \
    #            np.array( list[(6+(N_steps+1)*atom_2nd):(6+(N_steps+1)*atom_2nd+N_steps), 18], dtype=float )

#    plt.plot(E, dos1_s, label= str(atom_1st) + '_s')
 #   plt.plot(E, dos1_p_y, label= str(atom_1st) + '_p_y')
 #   plt.plot(E, dos1_p_z, label= str(atom_1st) + '_p_z')
 #   plt.plot(E, dos1_p_x, label= str(atom_1st) + '_p_x')
    plt.plot(E, dos1_d_xy, '--',  label= str(atom_1st) + '_d_xy')
    plt.plot(E, dos1_d_yz, '--', label= str(atom_1st) + '_d_yz')
    plt.plot(E, dos1_d_xz, '--', label= str(atom_1st) + '_d_xz')
    plt.plot(E, -dos1_d_z2, '--', label= str(atom_1st) + '_d_z2')
    plt.plot(E, -dos1_d_x2_y2, '--', label= str(atom_1st) + '_d_x2_y2')
    plt.plot(E, dos1_d_yz + dos1_d_xy + dos1_d_xz, 'k', label= str(atom_1st) + '_d_t2g')
    plt.plot(E, -dos1_d_z2 - dos1_d_x2_y2, 'k', label= str(atom_1st) + '_d_eg')
    plt.axhline(y=0)
    plt.axvline(x=0, ls='--')

#    plt.plot(E, -dos2_s, label= str(atom_2nd) + '_s')
 #   plt.plot(E, -dos2_p_y, label= str(atom_2nd) + '_p_y')
 #   plt.plot(E, -dos2_p_z, label= str(atom_2nd) + '_p_z')
 #   plt.plot(E, -dos2_p_x, label= str(atom_2nd) + '_p_x')
 #   plt.plot(E, -dos2_d_xy, label= str(atom_2nd) + '_d_xy')
 #   plt.plot(E, -dos2_d_yz, label= str(atom_2nd) + '_d_yz')
 #   plt.plot(E, -dos2_d_z2, label= str(atom_2nd) + '_d_z2')
 #   plt.plot(E, -dos2_d_xz, label= str(atom_2nd) + '_d_xz')
 #   plt.plot(E, -dos2_d_x2_y2, label= str(atom_2nd) + '_d_x2_y2')

    plt.legend(loc=0)
    plt.axis([axis_lim[0], axis_lim[1], axis_lim[2], axis_lim[3]])
    plt.xlabel('Energy (eV)')
    plt.ylabel('LDOS (State / atom / eV)')
    plt.savefig('LDOS.png')
 #   plt.show()


    # # spin-up
    # plt.plot(E, dos1_up_s, label= str(atom_1st) + '_s')
    # plt.plot(E, dos1_up_p, label= str(atom_1st) + '_p')
    # plt.plot(E, dos1_up_d, label= str(atom_1st) + '_d')
    # plt.plot(E, -dos2_up_s, label= str(atom_2nd) + '_s')
    # plt.plot(E, -dos2_up_p, label= str(atom_2nd) + '_p')
    # plt.plot(E, -dos2_up_d, label= str(atom_2nd) + '_d')
    # plt.legend(loc=0)
    # plt.axis([axis_lim[0], axis_lim[1], axis_lim[2]/2., axis_lim[3]/2.])
    # plt.xlabel('Energy (eV)')
    # plt.ylabel('LDOS (State / atom / eV)')
    # plt.savefig('LDOS-spin-up.png')
    # plt.close()

    # # spin-down
    # plt.plot(E, dos1_down_s, label= str(atom_1st) + '_s')
    # plt.plot(E, dos1_down_p, label= str(atom_1st) + '_p')
    # plt.plot(E, dos1_down_d, label= str(atom_1st) + '_d')
    # plt.plot(E, -dos2_down_s, label= str(atom_2nd) + '_s')
    # plt.plot(E, -dos2_down_p, label= str(atom_2nd) + '_p')
    # plt.plot(E, -dos2_down_d, label= str(atom_2nd) + '_d')
    # plt.legend(loc=0)
    # plt.axis([axis_lim[0], axis_lim[1], axis_lim[2]/2., axis_lim[3]/2.])
    # plt.xlabel('Energy (eV)')
    # plt.ylabel('LDOS (State / atom / eV)')
    # plt.savefig('LDOS-spin-down.png')
    # plt.close()

    # # spin-total
    # plt.plot(E, dos1_up_s + dos1_down_s, label= str(atom_1st) + '_s')
    # plt.plot(E, dos1_up_p + dos1_down_p, label= str(atom_1st) + '_p')
    # plt.plot(E, dos1_up_d + dos1_down_d, label= str(atom_1st) + '_d')
    # plt.plot(E, -dos2_up_s - dos2_down_s, label= str(atom_2nd) + '_s')
    # plt.plot(E, -dos2_up_p - dos2_down_p, label= str(atom_2nd) + '_p')
    # plt.plot(E, -dos2_up_d - dos2_down_d, label= str(atom_2nd) + '_d')
    # plt.legend(loc=0)
    # plt.axis(axis_lim)
    # plt.xlabel('Energy (eV)')
    # plt.ylabel('LDOS (State / atom / eV)')
    # plt.savefig('LDOS-spin-total.png')
    # plt.close()



#    table = np.column_stack((E, dos1_s, dos1_p, dos1_d, dos2_s, dos2_p, dos2_d))
#    np.savetxt('LDOS.txt', table, '%.6f', '\t')

elif is_spin and not is_m_decomposed:
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

elif not is_spin and is_m_decomposed:
    pass

elif not is_spin and not is_m_decomposed:
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
